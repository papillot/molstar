/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { BondType, ElementSymbol } from '../../../../mol-model/structure/model/types';
import { Unit } from '../../../../mol-model/structure/structure/unit';
import { type Structure } from '../../../../mol-model/structure/structure/structure';
import { type IntraUnitBonds } from '../../../../mol-model/structure/structure/unit/bonds/data';
import { type UnitIndex } from '../../../../mol-model/structure/structure/element/util';
import { eachBondedAtom, formalCharge, typeSymbol } from '../util';
import { type NumberArray } from '../../../../mol-util/type-helpers';
import { isDebugMode } from '../../../../mol-util/debug';
import { getElementIdx, isHydrogen } from '../../../../mol-model/structure/structure/unit/bonds/common';
import { Elements } from '../../../../mol-model/structure/model/properties/atomic/types';
import { OrderedSet } from '../../../../mol-data/int';

export interface State {
    structure: Structure
    unit: Unit.Atomic
    bonds: IntraUnitBonds
    start: UnitIndex
    end: UnitIndex
    /** UnitIndex per residue atom */
    unitIndices: UnitIndex[]
    /** element symbol per residue atom */
    el: ElementSymbol[]
    /** explicit hydrogen count per residue atom */
    hCount: Int8Array
    /** formal charge per residue atom */
    charge: Int8Array
    /** total covalent heavy degree (intra unit) per residue atom */
    degree: Int8Array
    /** intra-unit heavy neighbours (residue based index) per residue atom */
    heavyNeighbours: number[][]
    /** intra-unit all neighbours, explicit Hs included (residue based index) per residue atom */
    allNeighbours: number[][]
    /** geometry (from bond angles) per local atom, as `AtomGeometry` */
    geometry: Int8Array
    /** remaining open valence per local atom */
    open: Int8Array
    /** 1 if the atom belongs to a ring within the residue */
    inRing: Uint8Array
    /** `Ambiguity` bit-flags per local atom (geometry that can't be decided) */
    ambiguous: Uint8Array
    assignedBonds: Uint8Array,
    /** rings of size 5 or 6 that are fully contained in the residue, as local UnitIndex arrays */
    localRings5or6: number[][],
    /**
     * fused ring systems or isolated ring systems (ring components)
     * `atoms`: all atoms in the ring system.
     * `rings`: list of atoms per ring in the system.
     */
    localRingSystems: { atoms: number[], rings: number[][] }[],
}

export function State(structure: Structure, unit: Unit.Atomic, bonds: IntraUnitBonds, start: UnitIndex, end: UnitIndex): State {
    const { edgeProps, offset } = bonds;
    const { order, flags } = edgeProps;
    const n = end - start;

    const assignedBonds = new Uint8Array(order.length);
    const unitIndices = new Uint8Array(n) as unknown as UnitIndex[];
    const el = new Array<ElementSymbol>(n);
    const hCount = new Int8Array(n);
    const charge = new Int8Array(n);
    const degree = new Int8Array(n);
    const inRing = new Uint8Array(n);
    const heavyNeighbours: number[][] = [];
    const allNeighbours: number[][] = [];
    const indices = OrderedSet.ofRange(start, end);
    OrderedSet.forEach(indices, (u, i) => {
        unitIndices[i] = u;
        el[i] = typeSymbol(unit, u);
        const c = formalCharge(unit, u);
        charge[i] = Number.isFinite(c) ? c : 0;
        const heavyList: number[] = [];
        const allList: number[] = [];
        eachBondedAtom(structure, unit, u, (unitB, v) => {
            degree[i]++;
            if (unitB.id !== unit.id) return;
            const j = v - start;
            const ev = typeSymbol(unitB, v);
            if (isHydrogenElement(ev)) {
                hCount[i]++;
            } else {
                heavyList.push(j);
            }
            allList.push(j);
        });
        heavyNeighbours.push(heavyList);
        allNeighbours.push(allList);

        for (let j = offset[u]; j < offset[u + 1]; j++) {
            if (order[j] !== 1
                || (flags[j] & BondType.Flag.Computed) === 0
                || !BondType.isCovalent(flags[j])
            ) {
                assignedBonds[j] = 1;
            }
        }
    });

    const { all: allRings, elementRingIndices } = unit.rings;
    const seenRings = new Set<number>();
    // Rings of size 5 or 6 are special cases: collected
    const localRings5or6: number[][] = [];
    for (let u = start; u < end; u++) {
        const ringsForAtom = elementRingIndices.get(u);
        if (!ringsForAtom) continue;
        for (const ri of ringsForAtom) {
            if (seenRings.has(ri)) continue;
            seenRings.add(ri);
            const ring = allRings[ri];
            // ignore macrocycles spanning multiple residues
            if (ring[0] < start || ring[ring.length - 1] >= end) continue;
            const local = new Array<number>(ring.length);
            for (let i = 0; i < ring.length; i++) {
                const residueBasedIndex = ring[i] - start;
                local[i] = residueBasedIndex;
                inRing[residueBasedIndex] = 1;
            }
            if (ring.length < 5 || ring.length > 6) continue;
            const seq = ringTraversalSort(heavyNeighbours, local);
            if (seq) localRings5or6.push(seq);
        }
    }

    // Fused ring systems: group SSSR rings that share atoms (`unit.rings.ringComponents`). Keep only
    // components whose every member ring lies fully inside this residue and is size ≥5 — a component with
    // an out-of-residue or small (<5) member is skipped. Rings larger than 6 are kept and tested later
    // for aromaticity. A monocyclic ring is a 1-ring component.
    const localRingSystems: State['localRingSystems'] = [];
    for (const comp of unit.rings.ringComponents) {
        const atomsSet = new Set<number>();
        const compRings: number[][] = [];
        let ok = true;
        for (const ri of comp) {
            const ring = allRings[ri];
            if (ring[0] < start || ring[ring.length - 1] >= end) { ok = false; break; } // out of residue
            if (ring.length < 5) { ok = false; break; } // too small to be aromatic
            const local = new Array<number>(ring.length);
            for (let i = 0; i < ring.length; i++) { local[i] = ring[i] - start; atomsSet.add(local[i]); }
            compRings.push(local);
        }
        if (ok && atomsSet.size > 0) localRingSystems.push({ atoms: [...atomsSet], rings: compRings });
    }

    return {
        structure, unit, bonds, start, end, unitIndices,
        el, hCount, charge, degree, heavyNeighbours, allNeighbours, inRing,
        geometry: new Int8Array(n),
        open: new Int8Array(n),
        ambiguous: new Uint8Array(n),
        assignedBonds, localRings5or6, localRingSystems
    };
}

/**
 * order ring atoms indices by cyclic walk
 */
function ringTraversalSort(neighbours: number[][], ring: number[]): number[] | undefined {
    const len = ring.length;

    let current = ring[0], prev = -1;
    const seq = [current];
    for (let i = 1; i < len; i++) {
        for (const nb of neighbours[current]) {
            if (!ring.includes(nb) || prev === nb) continue;
            seq.push(nb);
            prev = current;
            current = nb;
            break;
        }
    }
    if (seq.length !== len) return undefined; // shoud never happen
    return seq;
}

/**
 * Cases where an atom's sp2/sp3 hybridization can't be decided from coordinates alone, so it
 * may take a double if a genuine sp2 neighbour and the bond geometry support it. Bit-flags so
 * more cases can be added later.
 */
export const enum Ambiguity {
    None = 0,
    /** degree-2 atom in a non-planar 5-ring: the ~108° angle can't tell sp2 from sp3 */
    Sp2OrSp3InRing5 = 1,
    /** degree-2 atom whose single angle sits just below the sp2 threshold (e.g. an aromatic-ring
     *  aldehyde C at 114.7°): too close to call sp2 vs sp3 from one angle */
    Sp2OrSp3Borderline = 2,
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();

/**
 * A bond is eligible for order perception if it is a single covalent bond whose order
 * is not authoritative. Such bonds are flagged `BondType.Flag.Computed`, which is set
 * both for distance-computed bonds and for PDB `CONECT` / `struct_conn` records that
 * give only basic connectivity without an explicit order (see `struct_conn.ts`). Bonds
 * with an explicit order (`chem_comp_bond`, struct_conn `doub`/`sing`, mol/mol2/sdf)
 * are not flagged `Computed` and are therefore left untouched.
 */
export function isPerceivable(flags: number, order: number) {
    return order === 1 && (flags & BondType.Flag.Computed) !== 0 && BondType.isCovalent(flags);
}

export function pos(state: State, local: number, out: Vec3) {
    const eI = state.unit.elements[state.start + local];
    return state.unit.conformation.position(eI, out);
}

export function distSq(state: State, aIdx: number, bIdx: number) {
    pos(state, aIdx, tmpVecA);
    pos(state, bIdx, tmpVecB);
    return Vec3.squaredDistance(tmpVecA, tmpVecB);
}



/** Set the order (and optionally OR-in flags) of the undirected bond u-v in both directed slots. */
export function setBond(bonds: IntraUnitBonds, u: UnitIndex, v: UnitIndex, order: number, addFlags: number, assignedBonds?: Uint8Array) {
    const ord = bonds.edgeProps.order as NumberArray;
    const flg = bonds.edgeProps.flags as NumberArray;
    let eidx = bonds.getDirectedEdgeIndex(u, v);
    if (eidx !== -1) {
        ord[eidx] = order;
        if (addFlags) flg[eidx] |= addFlags;
        if (assignedBonds) assignedBonds[eidx] = 1;
    }
    eidx = bonds.getDirectedEdgeIndex(v, u);
    if (eidx !== -1) {
        ord[eidx] = order;
        if (addFlags) flg[eidx] |= addFlags;
        if (assignedBonds) assignedBonds[eidx] = 1;
    }
    isDebugMode && console.log('Set bond ', u, '-', v, ' to order ', order, ' flags ', addFlags);
}

export function setAssignedBond(state: State, u: UnitIndex, v?: UnitIndex) {
    const { structure, unit, bonds } = state;
    if (v === void 0) {
        eachBondedAtom(structure, unit, u, (unitB, v) => {
            if (unitB.id !== unit.id) return;
            setAssignedBond(state, u, v);
        });
        return;
    }
    let eidx = bonds.getDirectedEdgeIndex(u, v);
    if (eidx !== -1) state.assignedBonds[eidx] = 1;
    eidx = bonds.getDirectedEdgeIndex(v, u);
    if (eidx !== -1) state.assignedBonds[eidx] = 1;
}

export function isAssignedBond(state: State, u: UnitIndex, v: UnitIndex) {
    const eidx = state.bonds.getEdgeIndex(u, v);
    if (eidx === -1) return true; // meaning we don't touch this bond
    return state.assignedBonds[eidx] === 1;
}

/** Whether atom `i` already has an incident multiple (double/triple) intra bond. */
export function hasMultipleBond(state: State, i: number) {
    const { unitIndices, heavyNeighbours } = state;
    const u = unitIndices[i];
    for (const j of heavyNeighbours[i]) {
        if (getOrder(state.bonds, u, unitIndices[j]) > 1) return true;
    }
    return false;
}

export function getOrder(bonds: IntraUnitBonds, u: UnitIndex, v: UnitIndex) {
    const eidx = bonds.getEdgeIndex(u, v);
    return eidx !== -1 ? bonds.edgeProps.order[eidx] : 0;
}

export function getFlags(bonds: IntraUnitBonds, u: UnitIndex, v: UnitIndex) {
    const eidx = bonds.getEdgeIndex(u, v);
    return eidx !== -1 ? bonds.edgeProps.flags[eidx] : BondType.Flag.None;
}

export function isHydrogenElement(el: ElementSymbol) {
    const elI = getElementIdx(el);
    return isHydrogen(elI);
}

/** Standard (Daylight-like) heavy valence used for open-valence bookkeeping. -1 = unknown. */
function defaultValence(el: ElementSymbol): number {
    switch (el) {
        case Elements.C: case Elements.SI: return 4;
        case Elements.N: case Elements.P: case Elements.B: return 3;
        case Elements.O: case Elements.S: case Elements.SE: return 2;
        case Elements.F: case Elements.CL: case Elements.BR: case Elements.I:
        case Elements.H: case Elements.D: case Elements.T: return 1;
        default: return -1;
    }
}

export function computeOpenValence(state: State) {
    const n = state.end - state.start;
    const { bonds, unitIndices } = state;
    for (let i = 0; i < n; i++) {
        let val = defaultValence(state.el[i]);
        if (val < 0) { state.open[i] = 0; continue; }
        val += state.charge[i];
        // sum current bond orders to heavy neighbours (intra + cross-residue counted as >=1)
        let used = state.hCount[i];
        const u = unitIndices[i];
        for (const j of state.heavyNeighbours[i]) {
            used += getOrder(bonds, u, unitIndices[j]);
        }
        state.open[i] = Math.max(0, val - used);
    }
}