from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, Spacer, Image, KeepTogether
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import inch, cm
from reportlab.lib.utils import ImageReader
import argparse
from draw_tree import get_example_tree

from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import HorizontalBarChart
from reportlab.lib.formatters import DecimalFormatter


class PdfWriter(object):
    def __init__(self):
        self.content = []
        self.small_space = KeepTogether(Spacer(0, 0.1 * inch))
        self.large_space = KeepTogether(Spacer(0, 0.3 * inch))
        self.styles = getSampleStyleSheet()

    @staticmethod
    def rgb_to_color_obj(rgb):
        decimal = tuple(float(x) / 255.0 for x in rgb)
        return colors.Color(red=decimal[0], green=decimal[1], blue=decimal[2])

    def draw_table(
            self,
            data,
            colwidths=None,
            cols_to_highlight=None,
            rows_to_highlight=None,
            border=True,
            grid=True,
            highlight_colour=(232, 237, 238)
    ):
        if cols_to_highlight is None:
            cols_to_highlight = list()
        if rows_to_highlight is None:
            rows_to_highlight = list()
        if colwidths is None:
            colwidths = '*'
        table = Table(data, hAlign='LEFT', colWidths=colwidths)
        style = []
        colour = self.rgb_to_color_obj(highlight_colour)
        if border:
            style.append(('BOX', (0, 0), (-1, -1), 0.25, colors.black))
        if grid:
            style.append(('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black))
        for col in cols_to_highlight:
            style.append(('BACKGROUND', (col, 0), (col, -1), colour))
            style.append(('FONTNAME', (col, 0), (col, -1), 'Helvetica-Bold'))
        for row in rows_to_highlight:
            style.append(('BACKGROUND', (0, row), (-1, row), colour))
            style.append(('FONTNAME', (0, row), (-1, row), 'Helvetica-Bold'))
        table.setStyle(TableStyle(style))
        self.content.append(table)
        self.content.append(self.large_space)

    def draw_title(self, text):
        self.content.append(Paragraph(text.upper(), self.styles['title']))
        self.content.append(self.large_space)

    def draw_heading(self, text, colour=(0, 94, 184), base_style='h4'):
        heading_style = ParagraphStyle(
            name='Modified_%s' % base_style,
            parent=self.styles[base_style],
            fontName='Helvetica-Bold',
            textColor=self.rgb_to_color_obj(colour)
        )
        heading = Paragraph(text.upper(), style=heading_style)
        self.content.append(heading)
        self.content.append(self.small_space)

    def draw_text(self, text):
        paragraph = Paragraph(text, style=self.styles['Normal'])
        self.content.append(paragraph)
        self.content.append(self.small_space)

    def plot_image(self, plot):
        img_reader = ImageReader(plot)
        iw, ih = img_reader.getSize()
        aspect = ih / float(iw)
        img = Image(plot, width=10 * cm, height=(10 * cm * aspect))
        self.content.append(img)

    def draw_bar_plot(self, data, names):
        d = Drawing(50, 100)
        bar = HorizontalBarChart()
        bar.data = data
        bar.x = 180
        bar.y = 0
        bar.categoryAxis.categoryNames = names
        bar.valueAxis.valueMin = 0
        bar.valueAxis.valueMax = 100
        bar.bars[0].fillColor = self.rgb_to_color_obj((65, 182, 230))
        bar.barLabels.angle = 0
        bar.barLabelFormat = DecimalFormatter(2, suffix='%')
        bar.barLabels.dx = 20
        d.add(bar, '')
        self.content.append(self.small_space)
        self.content.append(d)
        self.content.append(self.large_space)

    def compile_pdf(self, outfile):
        doc = SimpleDocTemplate(outfile)
        doc.build(self.content)


def parse_results(result_file, threshold):
    classified = 0
    results = []
    with open(result_file, 'r') as infile:
        lines = [l.strip().split('\t') for l in infile.readlines()]
    if lines:
        classified = int(lines[0][1])
        results_temp = lines[1:]
        for result in results_temp:
            if float(result[-1]) >= int(threshold):
                results.append(result)
        results.sort(key=lambda x: float(x[-1]), reverse=True)
    return classified, results


def classified_text(classified_reads, total_reads):
    """Converts read counts into sentence for adding to PDF document."""
    try:
        fraction = float(classified_reads) / float(total_reads)
        pct = '%.1f' % (fraction * 100)
    except ZeroDivisionError:
        pct = '0.0'
    text = 'Classified reads: {count} ({pct}%)'.format(count=classified_reads, pct=pct)
    return text


def parse_stats_file(filepath):
    with open(filepath) as infile:
        lines = [l.strip().split('\t') for l in infile.readlines()]
    stats = dict()
    for metric, score in lines:
        stats[metric] = score
    return stats


def write_results_to_pdf(outfile, seqmatch_output, centrifuge_output, stats_file, threshold=0):
    stats = parse_stats_file(stats_file)
    pdf_writer = PdfWriter()
    pdf_writer.draw_title('16s classification report')
    sample_data = [
        ['Run', stats.get('Run')],
        ['Sample ID', stats.get('Sample ID')],
        ['Sample Type', stats.get('Sample type')],
        ['Interval Time', ],
        ['Date', '']
    ]
    pdf_writer.draw_table(sample_data, cols_to_highlight=[0])

    analysed_reads = stats.get('Analysed reads')
    fmols = stats.get('Fmols in')
    if fmols:
        fmols = '%.2f' % fmols
    mean_read_length = stats.get('Mean read length')
    if mean_read_length:
        mean_read_length = '%.1f' % float(mean_read_length)
    mean_qscore = stats.get('Mean Q-score')
    if mean_qscore:
        mean_qscore = '%.2f' % float(mean_qscore)
    qual_data = [
        ['Total reads', stats.get('Total reads'), 'Fmols', fmols],
        ['Analysed reads', analysed_reads, 'Platform QC', stats.get('Platform QC')],
        ['Mean read length', mean_read_length, 'Mux Scan', stats.get('Mux Scan')],
        ['Mean Q-score', mean_qscore, 'Tapestation', stats.get('Tapestation')]
    ]
    pdf_writer.draw_heading('Quality summary')
    pdf_writer.draw_table(qual_data, cols_to_highlight=[0, 2])

    sqm_classified_reads, sqm_results = parse_results(seqmatch_output, threshold)
    cfg_classified_reads, cfg_results = parse_results(centrifuge_output, threshold)

    sqm_data = [['Species', 'Read Count', '%']]
    sqm_bar_data = []
    sqm_bar_names = []
    for result in sqm_results:
        sqm_data.append([result[0], result[1], result[2]])
        sqm_bar_data.append(float(result[2]))
        sqm_bar_names.append(result[0])

    pdf_writer.draw_heading('Seqmatch', base_style='h5')
    pdf_writer.draw_text(classified_text(sqm_classified_reads, analysed_reads))
    pdf_writer.draw_table(sqm_data, rows_to_highlight=[0], colwidths=[3.1*inch, 1.5*inch, 1.5*inch])
    sqm_bar_data.sort()
    pdf_writer.draw_bar_plot([sqm_bar_data], sqm_bar_names)

    cfg_data = [['Species', 'Rank', 'Read Count', '%']]
    cfg_tree_input_file = 'input_tree.tsv'
    with open(cfg_tree_input_file, 'w') as tree_input:
        tree_input.write('taxID\tsci_name\tread_prop\t#_reads')
        for result in cfg_results:
            cfg_data.append([result[0], result[2], result[3], result[4]])
            tree_input.write('\n{id}\t{name}\t{pct}\t{count}'.format(id=result[1], name=result[0], pct=result[4],
                                                                     count=result[3]))

    pdf_writer.draw_heading('Centrifuge/Supernatant', base_style='h5')
    pdf_writer.draw_text(classified_text(cfg_classified_reads, analysed_reads))
    pdf_writer.draw_table(cfg_data, rows_to_highlight=[0], colwidths=[2.8 * inch, 1.1 * inch, 1.1 * inch, 1.1 * inch])

    cfg_image = 'centrifuge.png'
    t, ts = get_example_tree(cfg_tree_input_file)
    t.render(cfg_image, w=1000, tree_style=ts)

    pdf_writer.plot_image(cfg_image)

    pdf_writer.compile_pdf(outfile)


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', required=True, action='store')
    parser.add_argument('--seqmatch_output', required=True, action='store')
    parser.add_argument('--centrifuge_output', required=True, action='store')
    parser.add_argument('--stats_file', required=True, action='store')
    parser.add_argument('--threshold', action='store', default=0)
    return parser.parse_args()


if __name__ == '__main__':
    args = argument_parser()
    write_results_to_pdf(args.output_file, args.seqmatch_output, args.centrifuge_output, args.stats_file,
                         int(args.threshold))
