/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Model, Unit } from '../../../../mol-model/structure';
import { BondType } from '../../../../mol-model/structure/model/types';
import { UnitIndex } from '../../../../mol-model/structure/structure/element/element';
import { IntraUnitBonds } from '../../../../mol-model/structure/structure/unit/bonds';
import { State, getFlags, getOrder, isPerceivable, setBond } from './common';

type BondingPattern = { a: string, b: string, order: number, flags: number }[];

const PerceivedCache = new WeakMap<Model, Map<string, BondingPattern>>();

function getModelCache(model: Model) {
    let c = PerceivedCache.get(model);
    if (!c) { c = new Map(); PerceivedCache.set(model, c); }
    return c;
}

function getCompHash(unit: Unit.Atomic, start: number, end: number, compId: string) {
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    const names: string[] = [];
    for (let i = start; i < end; i++) names.push(label_atom_id.value(unit.elements[i]));
    names.sort();
    return `${compId}|${names.join(',')}`;
}

function extractPattern(state: State): BondingPattern {
    const { unit, unitIndices, start, end, bonds } = state;
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    const pattern: BondingPattern = [];
    const n = end - start;
    for (let i = 0; i < n; i++) {
        const u = unitIndices[i];
        for (const j of state.heavyNeighbours[i]) {
            if (i >= j) continue;
            const v = unitIndices[j];
            const order = getOrder(bonds, u, v);
            const flags = getFlags(bonds, u, v);
            if (order > 1 || (flags & BondType.Flag.AromaticHuckel)) {
                pattern.push({
                    a: label_atom_id.value(unit.elements[u]),
                    b: label_atom_id.value(unit.elements[v]),
                    order,
                    flags: flags & (BondType.Flag.AromaticHuckel),
                });
            }
        }
    }
    return pattern;
}

export function applyCachedChemCompPattern(unit: Unit.Atomic, bonds: IntraUnitBonds, start: UnitIndex, end: UnitIndex, cachePrefix = ''): boolean {
    const model = unit.model;
    const cache = getModelCache(model);
    const { label_comp_id, label_atom_id } = model.atomicHierarchy.atoms;
    const compId = label_comp_id.value(unit.elements[start]);
    const compHash = cachePrefix + getCompHash(unit, start, end, compId);
        const pattern = cache.get(compHash);
    if (!pattern) return false;
    if (pattern.length === 0) return false;

    // TODO: handle altloc
    const nameToResidueBasedIndex = new Map<string, UnitIndex>();
    for (let i = start; i < end; i++) nameToResidueBasedIndex.set(label_atom_id.value(unit.elements[i]), i);
    for (const p of pattern) {
        const u = nameToResidueBasedIndex.get(p.a);
        const v = nameToResidueBasedIndex.get(p.b);
        if (u === undefined || v === undefined) continue;
        const flags = getFlags(bonds, u, v);
        if (!isPerceivable(flags, getOrder(bonds, u, v))) continue;
        setBond(bonds, u, v, p.order, BondType.Flag.Computed | p.flags);
    }
    return true;
}

export function cacheChemCompPattern(unit: Unit.Atomic, start: UnitIndex, end: UnitIndex, cachePrefix = '', state: State) {
    const model = unit.model;
    const cache = getModelCache(model);
    const { label_comp_id } = model.atomicHierarchy.atoms;
    const compId = label_comp_id.value(unit.elements[start]);
    const compHash = cachePrefix + getCompHash(unit, start, end, compId);
    if (cache.has(compHash)) return;
    const pattern = extractPattern(state);
    cache.set(compHash, pattern);
}
