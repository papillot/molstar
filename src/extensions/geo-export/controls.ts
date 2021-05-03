/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { PluginComponent } from '../../mol-plugin-state/component';
import { PluginContext } from '../../mol-plugin/context';
import { Task } from '../../mol-task';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateSelection } from '../../mol-state';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { SetUtils } from '../../mol-util/set';
import { ObjExporter } from './obj-exporter';
import { GlbExporter } from './glb-exporter';
import { StlExporter } from './stl-exporter';

export const GeometryParams = {
    format: PD.Select('obj', [['obj', 'OBJ + MTL'], ['glb', 'GLB'], ['stl', 'STL']])
};

export class GeometryControls extends PluginComponent {
    readonly behaviors = {
        params: this.ev.behavior<PD.Values<typeof GeometryParams>>(PD.getDefaultValues(GeometryParams))
    }

    private getFilename() {
        const models = this.plugin.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model)).map(s => s.obj!.data);
        const uniqueIds = new Set<string>();
        models.forEach(m => uniqueIds.add(m.entryId.toUpperCase()));
        const idString = SetUtils.toArray(uniqueIds).join('-');
        return `${idString || 'molstar-model'}`;
    }

    exportGeometries() {
        const task = Task.create('Export Geometries', async ctx => {
            try {
                const renderObjects = this.plugin.canvas3d?.getRenderObjects()!;
                const filename = this.getFilename();

                let renderObjectExporter: ObjExporter | GlbExporter | StlExporter;
                switch (this.behaviors.params.value.format) {
                    case 'obj':
                        renderObjectExporter = new ObjExporter(filename);
                        break;
                    case 'glb':
                        renderObjectExporter = new GlbExporter();
                        break;
                    case 'stl':
                        renderObjectExporter = new StlExporter();
                        break;
                    default: throw new Error('Unsupported format.');
                }

                for (let i = 0, il = renderObjects.length; i < il; ++i) {
                    await ctx.update({ message: `Exporting object ${i}/${il}` });
                    await renderObjectExporter.add(renderObjects[i], this.plugin.canvas3d?.webgl!, ctx);
                }

                const blob = await renderObjectExporter.getBlob(ctx);
                return {
                    blob,
                    filename: filename + '.' + renderObjectExporter.fileExtension
                };
            } catch (e) {
                this.plugin.log.error('' + e);
                throw e;
            }
        });

        return this.plugin.runTask(task, { useOverlay: true });
    }

    constructor(private plugin: PluginContext) {
        super();
    }
}