/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import T from '../tensor'
import { Mat4 } from '../linear-algebra-3d'

describe('tensor', () => {
    it('vector', () => {
        const V = T.Vector(3);
        const data = V.create();

        V.set(data, 0, 1);
        V.set(data, 1, 2);
        V.set(data, 2, 3);

        expect(data).toEqual(new Float64Array([1, 2, 3]));
        expect(V.get(data, 0)).toEqual(1);
        expect(V.get(data, 1)).toEqual(2);
        expect(V.get(data, 2)).toEqual(3);
    });

    it('matrix cm', () => {
        const M = T.ColumnMajorMatrix(3, 2, Int32Array);
        const data = M.create()

        // rows: [ [0, 1], [1, 2], [2, 3]  ]
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 2; j++) {
                M.set(data, i, j, i + j);
            }
        }

        expect(data).toEqual(new Int32Array([0, 1, 2, 1, 2, 3]));
    });

    it('matrix rm', () => {
        const M = T.RowMajorMatrix(3, 2, Int32Array);
        const data = M.create();

        // rows: [ [0, 1], [1, 2], [2, 3]  ]
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 2; j++) {
                M.set(data, i, j, i + j);
            }
        }
        expect(data).toEqual(new Int32Array([0, 1, 1, 2, 2, 3]));
    });

    it('mat4 equiv', () => {
        const M = T.ColumnMajorMatrix(4, 4);
        const data = M.create();
        const m = Mat4.zero();

        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                const v = (i + 1) * (j + 2);
                M.set(data, i, j, v);
                Mat4.setValue(m, i, j, v);
            }
        }

        for (let i = 0; i < 16; i++) expect(data[i]).toEqual(m[i]);
    });

    it('2d ij', () => {
        const M = T.Space([3, 4], [0, 1]);
        const data = M.create();
        const exp = new Float64Array(3 * 4)

        let o = 0;
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 4; j++) {
                M.set(data, i, j, o);
                exp[o] = o;
                o++;
            }
        }

        expect(data).toEqual(exp);
    });

    it('2d ji', () => {
        const M = T.Space([3, 4], [1, 0]);
        const data = M.create();
        const exp = new Float64Array(3 * 4)

        let o = 0;
        for (let j = 0; j < 4; j++) {
            for (let i = 0; i < 3; i++) {
                M.set(data, i, j, o);
                exp[o] = o;
                o++;
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d ijk', () => {
        const M = T.Space([3, 4, 5], [0, 1, 2]);
        const data = M.create();
        const exp = new Float64Array(3 * 4 * 5)

        let o = 0;
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 4; j++) {
                for (let k = 0; k < 5; k++) {
                    M.set(data, i, j, k, o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d ikj', () => {
        const M = T.Space([3, 3, 3], [0, 2, 1]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3)

        let o = 0;
        for (let i = 0; i < 3; i++) {
            for (let k = 0; k < 3; k++) {
                for (let j = 0; j < 3; j++) {
                    M.set(data, i, j, k, o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d jik', () => {
        const M = T.Space([3, 3, 3], [1, 0, 2]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3)

        let o = 0;
        for (let j = 0; j < 3; j++) {
            for (let i = 0; i < 3; i++) {
                for (let k = 0; k < 3; k++) {
                    M.set(data, i, j, k, o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });
    it('3d jki', () => {
        const M = T.Space([3, 3, 3], [1, 2, 0]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3)

        let o = 0;
        for (let j = 0; j < 3; j++) {
            for (let k = 0; k < 3; k++) {
                for (let i = 0; i < 3; i++) {
                    M.set(data, i, j, k, o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d kij', () => {
        const M = T.Space([3, 3, 3], [2, 0, 1]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3)

        let o = 0;
        for (let k = 0; k < 3; k++) {
            for (let i = 0; i < 3; i++) {
                for (let j = 0; j < 3; j++) {
                    M.set(data, i, j, k, o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d kji', () => {
        const M = T.Space([3, 3, 3], [2, 1, 0]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3)

        let o = 0;
        for (let k = 0; k < 3; k++) {
            for (let j = 0; j < 3; j++) {
                for (let i = 0; i < 3; i++) {
                    M.set(data, i, j, k, o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });
});