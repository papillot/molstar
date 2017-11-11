/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import SymmetryOperator from 'mol-math/geometry/symmetry-operator'
import AtomGroup from './atom/group'
import { Model } from '../model'

interface Unit extends SymmetryOperator.ArrayMapping {
    // Provides access to the underlying data.
    readonly model: Model,

    // Determines the operation applied to this unit.
    // The transform and and inverse are baked into the "getPosition" function
    readonly operator: SymmetryOperator,

    // The "full" atom group corresponding to this unit.
    readonly naturalGroup: AtomGroup,

    // Reference some commonly accessed things for faster access.
    readonly residueIndex: ArrayLike<number>,
    readonly chainIndex: ArrayLike<number>,
    readonly hierarchy: Model['hierarchy'],
    readonly conformation: Model['conformation']

    // TODO: add velocity?
}

namespace Unit {
    export function create(model: Model, operator: SymmetryOperator, naturalGroup: AtomGroup): Unit {
        const h = model.hierarchy;
        const { invariantPosition, position, x, y, z } = SymmetryOperator.createMapping(operator, model.conformation);

        return {
            model,
            operator,
            naturalGroup,
            residueIndex: h.residueSegments.segmentMap,
            chainIndex: h.chainSegments.segmentMap,
            hierarchy: model.hierarchy,
            conformation: model.conformation,
            invariantPosition,
            position,
            x, y, z
        };
    }

    export function withOperator(unit: Unit, operator: SymmetryOperator) {
        return create(unit.model, SymmetryOperator.compose(unit.operator, operator), unit.naturalGroup);
    }
}

export default Unit;