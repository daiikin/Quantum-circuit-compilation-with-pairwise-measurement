"""Render each diagram stage and the operational circuit without writing files."""

def render_stages(decoder, *, labels=False, show_checks=True, show_record_legend=False):
    import matplotlib.pyplot as plt
    import pyzx as zx
    from .flow import draw_flow
    from .circuit import draw_measurement_circuit

    compilation = decoder.compilation
    measurement = compilation.measurement
    signed = ('-' if measurement.sign < 0 else '+') + measurement.pauli
    figures = {}
    try:
        def zx_figure(key, graph, title):
            fig = zx.draw_matplotlib(graph, labels=labels)
            figures[key] = fig
            fig.suptitle(title)

        zx_figure('initial', measurement.original, f'Initial ZX branch: {signed}')
        if measurement.rewritten is not None:
            zx_figure('direct_rewrite', measurement.rewritten, 'Optional direct R4/R5 preview')
        zx_figure('prepared', compilation.prepared.graph, 'After MC-flow preparation')
        # Flow rendering lays out a copy; it does not mutate the compiled graph.
        figures['reduced'] = draw_flow(compilation.reduced.graph, compilation.reduced.paths,
                                       labels=labels, title='Reduced ZX diagram with MC-flow layout')
        figures['circuit'] = draw_measurement_circuit(
            decoder, title=f'Operational measurement of {signed}',
            show_checks=show_checks, show_record_legend=show_record_legend)
        return figures
    except Exception:
        for fig in figures.values():
            plt.close(fig)
        raise
