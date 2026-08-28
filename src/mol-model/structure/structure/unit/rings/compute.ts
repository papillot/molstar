/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Segmentation, SortedArray } from '../../../../../mol-data/int';
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import { BondType } from '../../../model/types';
import { StructureElement } from '../../../structure';
import { Unit } from '../../unit';
import { IntraUnitBonds } from '../bonds/data';
import { sortArray } from '../../../../../mol-data/util';
import { Column } from '../../../../../mol-data/db';
import { arraySetAdd, arraySetRemove } from '../../../../../mol-util/array';

export function computeRings(unit: Unit.Atomic, bonds: IntraUnitBonds = unit.bonds) {
    const size = largestResidue(unit);
    const state = State(unit, size, bonds);

    const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
    while (residuesIt.hasNext) {
        const seg = residuesIt.move();
        processResidue(state, seg.start, seg.end);
    }

    return state.rings;
}

const enum Constants {
    MaxDepth = 5
}

interface State {
    startVertex: number,
    endVertex: number,
    count: number,
    isRingAtom: Int32Array,
    marked: Int32Array,
    queue: Int32Array,
    color: Int32Array,
    pred: Int32Array,
    depth: Int32Array,

    left: Int32Array,
    right: Int32Array,

    currentColor: number,
    currentAltLoc: string,
    hasAltLoc: boolean,

    rings: SortedArray<StructureElement.UnitIndex>[],
    currentRings: SortedArray<StructureElement.UnitIndex>[],
    bonds: IntraUnitBonds,
    unit: Unit.Atomic,
    altLoc: Column<string>
}

function State(unit: Unit.Atomic, capacity: number, bonds: IntraUnitBonds): State {
    return {
        startVertex: 0,
        endVertex: 0,
        count: 0,
        isRingAtom: new Int32Array(capacity),
        marked: new Int32Array(capacity),
        queue: new Int32Array(capacity),
        pred: new Int32Array(capacity),
        depth: new Int32Array(capacity),
        left: new Int32Array(Constants.MaxDepth),
        right: new Int32Array(Constants.MaxDepth),
        color: new Int32Array(capacity),
        currentColor: 0,
        currentAltLoc: '',
        hasAltLoc: false,
        rings: [],
        currentRings: [],
        unit,
        bonds,
        altLoc: unit.model.atomicHierarchy.atoms.label_alt_id
    };
}

function resetState(state: State) {
    state.count = state.endVertex - state.startVertex;
    const { isRingAtom, pred, color, depth, marked } = state;
    for (let i = 0; i < state.count; i++) {
        isRingAtom[i] = 0;
        pred[i] = -1;
        marked[i] = -1;
        color[i] = 0;
        depth[i] = 0;
    }
    state.currentColor = 0;
    state.currentAltLoc = '';
    state.hasAltLoc = false;
}

function resetDepth(state: State) {
    const { depth } = state;
    for (let i = 0; i < state.count; i++) {
        depth[i] = state.count + 1;
    }
}

function largestResidue(unit: Unit.Atomic) {
    const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
    let size = 0;
    while (residuesIt.hasNext) {
        const seg = residuesIt.move();
        size = Math.max(size, seg.end - seg.start);
    }
    return size;
}

function isStartIndex(state: State, i: number) {
    const bondOffset = state.bonds.offset;
    const a = state.startVertex + i;
    const bStart = bondOffset[a], bEnd = bondOffset[a + 1];
    const bondCount = bEnd - bStart;
    if (bondCount <= 1 || (state.isRingAtom[i] && bondCount === 2)) return false;
    return true;
}

function processResidue(state: State, start: number, end: number) {
    state.startVertex = start;
    state.endVertex = end;

    // no two atom rings
    if (state.endVertex - state.startVertex < 3) return;

    state.currentRings = [];

    const { elements } = state.unit;
    const altLocs: string[] = [];
    for (let i = state.startVertex; i < state.endVertex; i++) {
        const altLoc = state.altLoc.value(elements[i]);
        arraySetAdd(altLocs, altLoc);
    }
    arraySetRemove(altLocs, '');

    let mark = 1;
    if (altLocs.length === 0) {
        resetState(state);
        for (let i = 0; i < state.count; i++) {
            if (!isStartIndex(state, i)) continue;
            resetDepth(state);
            mark = findRings(state, i, mark);
        }
        recoverLargeCycles(state);
    } else {
        for (let aI = 0; aI < altLocs.length; aI++) {
            resetState(state);
            state.hasAltLoc = true;
            state.currentAltLoc = altLocs[aI];
            for (let i = 0; i < state.count; i++) {
                if (!isStartIndex(state, i)) continue;
                const altLoc = state.altLoc.value(elements[state.startVertex + i]);
                if (altLoc && altLoc !== state.currentAltLoc) {
                    continue;
                }
                resetDepth(state);
                mark = findRings(state, i, mark);
            }
        }
    }

    for (let i = 0, _i = state.currentRings.length; i < _i; i++) {
        state.rings.push(state.currentRings[i]);
    }
}

function addRing(state: State, a: number, b: number, isRingAtom: Int32Array) {
    // only "monotonous" rings
    if (b < a) {
        return false;
    }

    const { pred, color, left, right } = state;
    const nc = ++state.currentColor;

    let current = a;

    for (let t = 0; t < Constants.MaxDepth; t++) {
        color[current] = nc;
        current = pred[current];
        if (current < 0) break;
    }

    let leftOffset = 0, rightOffset = 0;

    let found = false, target = 0;
    current = b;
    for (let t = 0; t < Constants.MaxDepth; t++) {
        if (color[current] === nc) {
            target = current;
            found = true;
            break;
        }
        right[rightOffset++] = current;
        current = pred[current];
        if (current < 0) break;
    }
    if (!found) {
        return false;
    }

    current = a;
    for (let t = 0; t < Constants.MaxDepth; t++) {
        left[leftOffset++] = current;
        if (target === current) break;
        current = pred[current];
        if (current < 0) break;
    }

    const len = leftOffset + rightOffset;
    // rings must have at least three elements
    if (len < 3) {
        return false;
    }

    const ring = new Int32Array(len);
    let ringOffset = 0;
    for (let t = 0; t < leftOffset; t++) {
        ring[ringOffset++] = state.startVertex + left[t];
        isRingAtom[left[t]] = 1;
    }
    for (let t = rightOffset - 1; t >= 0; t--) {
        ring[ringOffset++] = state.startVertex + right[t];
        isRingAtom[right[t]] = 1;
    }

    sortArray(ring);

    // Check if the ring is unique and another one is not it's subset
    for (let rI = 0, _rI = state.currentRings.length; rI < _rI; rI++) {
        const r = state.currentRings[rI];

        if (ring.length === r.length) {
            if (SortedArray.areEqual(ring as any, r)) return false;
        } else if (ring.length > r.length) {
            if (SortedArray.isSubset(ring as any, r)) return false;
        }
    }

    state.currentRings.push(SortedArray.ofSortedArray(ring));

    return true;
}

function findRings(state: State, from: number, mark: number) {
    const { bonds, startVertex, endVertex, isRingAtom, marked, queue, pred, depth } = state;
    const { elements } = state.unit;
    const { b: neighbor, edgeProps: { flags: bondFlags }, offset } = bonds;
    marked[from] = mark;
    depth[from] = 0;
    queue[0] = from;
    let head = 0, size = 1;

    while (head < size) {
        const top = queue[head++];
        const d = depth[top];
        const a = startVertex + top;
        const start = offset[a], end = offset[a + 1];

        for (let i = start; i < end; i++) {
            const b = neighbor[i];
            if (b < startVertex || b >= endVertex || !BondType.isCovalent(bondFlags[i])) continue;

            if (state.hasAltLoc) {
                const altLoc = state.altLoc.value(elements[b]);
                if (altLoc && state.currentAltLoc !== altLoc) {
                    continue;
                }
            }

            const other = b - startVertex;

            if (marked[other] === mark) {
                if (pred[other] !== top && pred[top] !== other) {
                    if (addRing(state, top, other, isRingAtom)) {
                        return mark + 1;
                    }
                }
                continue;
            }

            const newDepth = Math.min(depth[other], d + 1);
            if (newDepth > Constants.MaxDepth) continue;

            depth[other] = newDepth;
            marked[other] = mark;
            queue[size++] = other;
            pred[other] = top;
        }
    }
    return mark + 1;
}

/** Residue-local covalent adjacency (both endpoints in `[startVertex, endVertex)`); indices are `atom - startVertex`. */
function buildResidueAdjacency(state: State, count: number): number[][] {
    const { bonds, startVertex, endVertex } = state;
    const { b: neighbor, edgeProps: { flags }, offset } = bonds;
    const adj: number[][] = new Array(count);
    for (let i = 0; i < count; i++) adj[i] = [];
    for (let i = 0; i < count; i++) {
        const a = startVertex + i;
        for (let j = offset[a], jl = offset[a + 1]; j < jl; j++) {
            const nb = neighbor[j];
            if (nb <= a || nb < startVertex || nb >= endVertex) continue; // each undirected edge once
            if (!BondType.isCovalent(flags[j])) continue;
            const other = nb - startVertex;
            adj[i].push(other);
            adj[other].push(i);
        }
    }
    return adj;
}

/**
 * Recover cycles too large for the bounded BFS (which cannot emit rings > 2 * MaxDepth atoms).
 * Euler gate: in each connected component of the residue covalent graph the number of independent
 * rings equals E - V + 1 (cyclomatic number); if the BFS found fewer, larger cycles are provably
 * missing and are recovered with a Horton minimum cycle basis. Runs only on a cyclomatic mismatch,
 * so ~every residue (proteins, sugars, ordinary fused systems) is untouched.
 */
function recoverLargeCycles(state: State) {
    const count = state.endVertex - state.startVertex;
    // an undetected cycle is necessarily larger than 2 * MaxDepth atoms, so a component must have at
    // least that many atoms; skip smaller residues cheaply before doing any graph analysis.
    if (count <= 2 * Constants.MaxDepth) return;

    const adj = buildResidueAdjacency(state, count);

    // connected components of the covalent graph (iterative flood fill)
    const comp = new Int32Array(count).fill(-1);
    let nComp = 0;
    const stack: number[] = [];
    for (let s = 0; s < count; s++) {
        if (comp[s] !== -1 || adj[s].length === 0) continue;
        comp[s] = nComp;
        stack.length = 0; stack.push(s);
        while (stack.length) {
            const v = stack.pop()!;
            for (const w of adj[v]) if (comp[w] === -1) { comp[w] = nComp; stack.push(w); }
        }
        nComp++;
    }
    if (nComp === 0) return;

    // per component: vertices, undirected edge count, and rings already found by the BFS
    const compVerts: number[][] = [];
    const compEdges = new Int32Array(nComp);
    const compFound: SortedArray<StructureElement.UnitIndex>[][] = [];
    for (let c = 0; c < nComp; c++) { compVerts.push([]); compFound.push([]); }
    for (let v = 0; v < count; v++) {
        if (comp[v] === -1) continue;
        compVerts[comp[v]].push(v);
        compEdges[comp[v]] += adj[v].length;
    }
    for (let c = 0; c < nComp; c++) compEdges[c] >>= 1; // each undirected edge counted from both ends
    for (const ring of state.currentRings) {
        const c = comp[ring[0] - state.startVertex];
        if (c >= 0) compFound[c].push(ring);
    }

    for (let c = 0; c < nComp; c++) {
        const mu = compEdges[c] - compVerts[c].length + 1; // cyclomatic number of the component
        if (compFound[c].length >= mu) continue; // BFS already found every independent ring
        hortonRecover(state, adj, compVerts[c], compFound[c], mu);
    }
}

/** Set bit `i` in a Uint32Array bitset. */
function setBit(vec: Uint32Array, i: number) { vec[i >> 5] |= (1 << (i & 31)); }
/** Highest set bit index, or -1 if the bitset is zero. */
function highestBit(vec: Uint32Array): number {
    for (let w = vec.length - 1; w >= 0; w--) {
        if (vec[w] !== 0) return w * 32 + (31 - Math.clz32(vec[w]));
    }
    return -1;
}

/**
 * Horton minimum cycle basis on one connected component. Seeds a GF(2) edge-incidence basis with the
 * rings the BFS already found, then adds the shortest Horton candidate cycles that increase the rank
 * until it reaches `mu` (the cyclomatic number). The candidates added beyond the seeds are the missing
 * minimal cycles (e.g. the porphyrin 16-macrocycle); they are appended to `state.currentRings`.
 */
function hortonRecover(state: State, adj: number[][], verts: number[], foundRings: SortedArray<StructureElement.UnitIndex>[], mu: number) {
    // reduce to the 2-core: iteratively drop degree <= 1 vertices (side chains carry no cycles)
    const deg = new Map<number, number>();
    for (const v of verts) deg.set(v, adj[v].length);
    const removed = new Set<number>();
    let changed = true;
    while (changed) {
        changed = false;
        for (const v of verts) {
            if (removed.has(v) || deg.get(v)! > 1) continue;
            removed.add(v); changed = true;
            for (const w of adj[v]) if (!removed.has(w)) deg.set(w, deg.get(w)! - 1);
        }
    }
    const core: number[] = [];
    for (const v of verts) if (!removed.has(v)) core.push(v);
    const k = core.length;
    if (k === 0) return;

    // compact core indices 0..k-1 and their adjacency
    const idx = new Map<number, number>();
    for (let i = 0; i < k; i++) idx.set(core[i], i);
    const cadj: number[][] = new Array(k);
    for (let i = 0; i < k; i++) cadj[i] = [];
    for (let i = 0; i < k; i++) {
        for (const w of adj[core[i]]) {
            const wi = idx.get(w);
            if (wi !== undefined) cadj[i].push(wi);
        }
    }

    // undirected edge indexing
    const edgeId = new Map<number, number>();
    let nEdges = 0;
    for (let a = 0; a < k; a++) for (const b of cadj[a]) if (a < b) edgeId.set(a * k + b, nEdges++);
    const eid = (a: number, b: number) => edgeId.get(a < b ? a * k + b : b * k + a)!;
    const W = (nEdges + 31) >> 5;

    // all-pairs shortest paths on the core (BFS from each vertex)
    const dist: Int32Array[] = new Array(k);
    const bfsQueue = new Int32Array(k);
    for (let s = 0; s < k; s++) {
        const d = new Int32Array(k).fill(-1);
        let head = 0, tail = 0; bfsQueue[tail++] = s; d[s] = 0;
        while (head < tail) {
            const u = bfsQueue[head++];
            for (const w of cadj[u]) if (d[w] < 0) { d[w] = d[u] + 1; bfsQueue[tail++] = w; }
        }
        dist[s] = d;
    }
    // shortest path s -> t as a vertex list [t, ..., s], or null if unreachable
    const path = (s: number, t: number): number[] | null => {
        const ds = dist[s];
        if (ds[t] < 0) return null;
        const p = [t];
        let cur = t;
        while (cur !== s) {
            let next = -1;
            for (const w of cadj[cur]) if (ds[w] === ds[cur] - 1) { next = w; break; }
            if (next < 0) return null;
            p.push(next); cur = next;
        }
        return p;
    };

    // GF(2) basis keyed by leading bit; returns true if the vector was independent (and added)
    const basisByLead = new Map<number, Uint32Array>();
    const addToBasis = (vec: Uint32Array): boolean => {
        let lead = highestBit(vec);
        while (lead >= 0) {
            const existing = basisByLead.get(lead);
            if (!existing) { basisByLead.set(lead, vec); return true; }
            for (let w = 0; w < W; w++) vec[w] ^= existing[w];
            lead = highestBit(vec);
        }
        return false;
    };

    // seed with the rings already found (chordless minimal rings → induced edges are the ring edges)
    for (const ring of foundRings) {
        const rverts: number[] = [];
        let ok = true;
        for (let i = 0; i < ring.length; i++) {
            const ci = idx.get(ring[i] - state.startVertex);
            if (ci === undefined) { ok = false; break; }
            rverts.push(ci);
        }
        if (!ok) continue;
        const rset = new Set(rverts);
        const vec = new Uint32Array(W);
        for (const a of rverts) for (const b of cadj[a]) if (a < b && rset.has(b)) setBit(vec, eid(a, b));
        addToBasis(vec);
    }
    if (basisByLead.size >= mu) return;

    // Horton candidates: for each vertex v and edge (x,y), if paths v->x and v->y meet only at v,
    // the cycle v..x-y..v has length d(v,x)+d(v,y)+1. Collect, sort ascending, add the independent ones.
    type Candidate = { len: number, verts: number[], vec: Uint32Array };
    const candidates: Candidate[] = [];
    for (let v = 0; v < k; v++) {
        for (let x = 0; x < k; x++) {
            for (const y of cadj[x]) {
                if (x >= y) continue; // each edge once
                if (v === x || v === y) continue; // degenerate: cycle would be the edge itself
                const px = path(v, x), py = path(v, y);
                if (!px || !py) continue;
                const seen = new Set(px);
                let disjoint = true;
                for (const p of py) if (p !== v && seen.has(p)) { disjoint = false; break; }
                if (!disjoint) continue;
                const cverts = new Set<number>(px);
                for (const p of py) cverts.add(p);
                if (cverts.size < 3) continue; // not a real ring
                const vec = new Uint32Array(W);
                for (let i = 0; i + 1 < px.length; i++) setBit(vec, eid(px[i], px[i + 1]));
                for (let i = 0; i + 1 < py.length; i++) setBit(vec, eid(py[i], py[i + 1]));
                setBit(vec, eid(x, y));
                candidates.push({ len: (px.length - 1) + (py.length - 1) + 1, verts: [...cverts], vec });
            }
        }
    }
    candidates.sort((a, b) => a.len - b.len);

    for (const cand of candidates) {
        if (basisByLead.size >= mu) break;
        if (!addToBasis(cand.vec)) continue; // dependent on what we already have
        // emit the recovered ring in unit-index space
        const ring = new Int32Array(cand.verts.length);
        for (let i = 0; i < cand.verts.length; i++) ring[i] = state.startVertex + core[cand.verts[i]];
        sortArray(ring);
        appendRecoveredRing(state, ring);
    }
}

/** Append a recovered ring to `currentRings`, applying the same uniqueness/subset filter as `addRing`. */
function appendRecoveredRing(state: State, ring: Int32Array) {
    for (let rI = 0, _rI = state.currentRings.length; rI < _rI; rI++) {
        const r = state.currentRings[rI];
        if (ring.length === r.length) {
            if (SortedArray.areEqual(ring as any, r)) return;
        } else if (ring.length > r.length) {
            if (SortedArray.isSubset(ring as any, r)) return;
        }
    }
    state.currentRings.push(SortedArray.ofSortedArray(ring));
}

export function getFingerprint(elements: string[]) {
    const len = elements.length;
    const reversed: string[] = new Array(len);

    for (let i = 0; i < len; i++) reversed[i] = elements[len - i - 1];

    const rotNormal = getMinimalRotation(elements);
    const rotReversed = getMinimalRotation(reversed);

    let isNormalSmaller = false;

    for (let i = 0; i < len; i++) {
        const u = elements[(i + rotNormal) % len], v = reversed[(i + rotReversed) % len];
        if (u !== v) {
            isNormalSmaller = u < v;
            break;
        }
    }

    if (isNormalSmaller) return buildFinderprint(elements, rotNormal);
    return buildFinderprint(reversed, rotReversed);
}

function getMinimalRotation(elements: string[]) {
    // adapted from http://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation

    const len = elements.length;
    const f = new Int32Array(len * 2);
    for (let i = 0; i < f.length; i++) f[i] = -1;

    let u = '', v = '', k = 0;

    for (let j = 1; j < f.length; j++) {
        let i = f[j - k - 1];
        while (i !== -1) {
            u = elements[j % len]; v = elements[(k + i + 1) % len];
            if (u === v) break;
            if (u < v) k = j - i - 1;
            i = f[i];
        }

        if (i === -1) {
            u = elements[j % len]; v = elements[(k + i + 1) % len];
            if (u !== v) {
                if (u < v) k = j;
                f[j - k] = -1;
            } else f[j - k] = i + 1;
        } else f[j - k] = i + 1;
    }

    return k;
}

function buildFinderprint(elements: string[], offset: number) {
    const len = elements.length;
    const ret: string[] = [];
    let i;
    for (i = 0; i < len - 1; i++) {
        ret.push(elements[(i + offset) % len]);
        ret.push('-');
    }
    ret.push(elements[(i + offset) % len]);
    return ret.join('');
}

type RingIndex = import('../rings').UnitRings.Index
type RingComponentIndex = import('../rings').UnitRings.ComponentIndex

export function createIndex(rings: ArrayLike<SortedArray<StructureElement.UnitIndex>>, aromaticRings: ReadonlyArray<RingIndex>) {
    const elementRingIndices: Map<StructureElement.UnitIndex, RingIndex[]> = new Map();
    const elementAromaticRingIndices: Map<StructureElement.UnitIndex, RingIndex[]> = new Map();

    // for each ring atom, assign all rings that it is present in
    for (let rI = 0 as RingIndex, _rI = rings.length; rI < _rI; rI++) {
        const r = rings[rI];
        for (let i = 0, _i = r.length; i < _i; i++) {
            const e = r[i];
            if (elementRingIndices.has(e)) elementRingIndices.get(e)!.push(rI);
            else elementRingIndices.set(e, [rI]);
        }
    }

    // for each ring atom, assign all aromatic rings that it is present in
    for (let aI = 0, _aI = aromaticRings.length; aI < _aI; aI++) {
        const rI = aromaticRings[aI];
        const r = rings[rI];
        for (let i = 0, _i = r.length; i < _i; i++) {
            const e = r[i];
            if (elementAromaticRingIndices.has(e)) elementAromaticRingIndices.get(e)!.push(rI);
            else elementAromaticRingIndices.set(e, [rI]);
        }
    }

    // create a graph where vertices are rings, edge if two rings share at least one atom
    const graph = new IntAdjacencyGraph.UniqueEdgeBuilder(rings.length);
    for (let rI = 0 as RingIndex, _rI = rings.length; rI < _rI; rI++) {
        const r = rings[rI];

        for (let i = 0, _i = r.length; i < _i; i++) {
            const e = r[i];

            const containedRings = elementRingIndices.get(e)!;

            if (containedRings.length === 1) continue;

            for (let j = 0, _j = containedRings.length; j < _j; j++) {
                const rJ = containedRings[j];
                if (rI >= rJ) continue;
                graph.addEdge(rI, rJ);
            }
        }
    }

    const components = IntAdjacencyGraph.connectedComponents(graph.getGraph());

    const ringComponentIndex = components.componentIndex as any as RingComponentIndex[];
    const ringComponents: RingIndex[][] = [];
    for (let i = 0; i < components.componentCount; i++) ringComponents[i] = [];

    for (let rI = 0 as RingIndex, _rI = rings.length; rI < _rI; rI++) {
        ringComponents[ringComponentIndex[rI]].push(rI);
    }

    return { elementRingIndices, elementAromaticRingIndices, ringComponentIndex, ringComponents };
}