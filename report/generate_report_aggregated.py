import argparse
import logging
import pandas as pd
import yaml
import re
from pathlib import Path
import numpy as np
from rich.logging import RichHandler
from typing import Dict, Any, Optional, Sequence

#############
# Jinja env #
#############
from jinja2 import Environment, FileSystemLoader

env = Environment(loader=FileSystemLoader(Path(__file__).absolute().parent / "report_template"))
template = env.get_template("report_aggregated.html")


def generate_plots(output_path: Path, input_path: Path, config_data: Dict[str, Any]) -> None:
    """
    Runs snakemake.
    :param output_path: Directory where data and final report is placed
    :param input_path: Directory where all output.tsv files of the classifieres reside.
    :param config_data: Arguments of this function which are saved to a file for later reference.
    :return:None
    """
    import analysis_scripts.analysis_aggregated
    # Create output directory
    output_path_data = output_path / 'data'
    if not output_path_data.exists():
        output_path_data.mkdir(parents=True)

    # Create config
    config_path = output_path / 'config.yml'
    with config_path.open('w') as handle:
        yaml.dump(config_data, handle)
    logger.info(f"Config created: {config_path}")

    # Execute
    analysis_scripts.analysis_aggregated.main(input_path, output_path_data)


def generate_report(output_path: Path) -> None:
    template_vars = {"title": 'Aggregated report'}
    output_path_data = output_path / 'data'

    with open(output_path_data / 'used_datasets.txt', 'r') as f:
        used_datasets = f.read().rstrip().replace('\n', ', ')

    template_vars['used_datasets'] = used_datasets

    images_names = ['boxplot_metrics',
                    'cutoff_metrics_plot',
                    'cutoff_metrics_plot_Q',
                    'PR_errorbars'
                    ]
    table_names = ['metric_table']

    def color_cells(s):
        bins = np.linspace(0.1, 1, 9)
        colors = ['rgb(255,247,251)', 'rgb(236,226,240)', 'rgb(208,209,230)', 'rgb(166,189,219)', 'rgb(103,169,207)', 'rgb(54,144,192)',
                  'rgb(2,129,138)', 'rgb(1,108,89)', 'rgb(1,70,54)']
        bin_number = np.digitize(s, bins)
        if s == 0 or np.isnan(s):
            return None
        elif s >= 1:
            return f'background-color: rgb(0, 37, 27); color: white;'
        else:
            if bin_number >=7:
                return f'background-color: {colors[bin_number]}; color: white;'
            else:
                return f'background-color: {colors[bin_number]}'

    # images = {}
    for organism in ['genus', 'species']:
        # images[organism] = {}
        for plot in images_names:
            with open(output_path_data / organism / (plot + '.svg'), 'r') as svg:
                image = svg.read()
                image = image[image.find('<svg'):]
                template_vars[plot + '_' + organism] = image
        for table in table_names:
            table_path = output_path_data / organism / (table + '.tsv')
            df = (pd.read_table(table_path, index_col=[0, 1], sep='\t').style.applymap(color_cells)
                  .set_table_attributes('border="3" class="table-condensed"')  # class is in jinja template
                  .to_html(float_format='%.5g'))
            template_vars['table_style' + '_' + organism] = re.search('<style type="text/css">\n(.*)</style>',
                                                                      df,
                                                                      flags=re.DOTALL).group(1)
            template_vars[table + '_' + organism] = re.search('<table.*>\n(.*)</table>',
                                                              df,
                                                              flags=re.DOTALL).group()

    html_out = template.render(template_vars)
    with open(output_path / 'generated_report.html', 'w+') as f:
        f.write(html_out)


class Report(object):
    """
    This class contains is used to generate the report .
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the a report instance.
        """
        self._args = Report._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', metavar='', help="Output from SnakeMake flow", required=True)
        parser.add_argument('-o', '--output', metavar='', help='Location of report/figures', required=True)

        return parser.parse_args(args)

    def __get_config_data(self) -> Dict[str, Any]:
        """
        Get the snakemake configuration data.
        :return: Configuration data
        """
        return {
            'output': self._args.output,
            'input': self._args.input,
        }

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """

        self.__get_config_data()
        generate_plots(Path(self._args.output), Path(self._args.input), self.__get_config_data())
        logger.info(f"Images created")

        generate_report(Path(self._args.output))


if __name__ == '__main__':
    # create logger
    logger = logging.getLogger('classification_agg')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(RichHandler(show_path=False, level=logging.INFO,
                                  ))

    report = Report()
    report.run()
