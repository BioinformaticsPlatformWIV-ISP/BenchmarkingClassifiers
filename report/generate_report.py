import argparse
import logging
import sys
import pandas as pd
from rich.logging import RichHandler
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, Sequence

#############
# Jinja env #
#############
from jinja2 import Environment, FileSystemLoader

REPORT_TEMPLATE = Path(__file__).parent / 'report_template'

env = Environment(loader=FileSystemLoader(REPORT_TEMPLATE))
template = env.get_template("report_css.html")


def generate_plots(output_path: Path, input_path: Path, classifiers, config_data: Dict[str, Any]) -> None:
    """
    Runs snakemake.
    :param output_path: Directory where data and final report is placed
    :param input_path: Directory where all output.tsv files of the classifieres reside.
    :param config_data: Arguments of this function which are saved to a file for later reference.
    :param classifiers: For which classifiers should a report be generated.
    :return:None
    """
    # Before executing any command, check if structure of input folder is correct
    if not all([(input_path / 'genus').exists(), (input_path / 'species').exists()]):
        logger.error("One of 'genus' or 'species' subfolders does not exist. Exiting...")
        sys.exit()

    # Before executing any command, check if a ground truth file is present in the input folder
    if not (input_path.parents[0] / 'ground_truth_updated.yml').exists():
        logger.error(f'One of "ground_truth_updated.yml" or "ground_truth_updated.yml" not found in {input_path.resolve()}. Exiting..."')
        exit()

    import analysis_scripts.analysis
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
    analysis_scripts.analysis.main(input_path, output_path_data, classifiers)


def generate_report(output_path: Path, title) -> None:
    template_vars = {"title": title}
    output_path_data = output_path / 'data'

    images_names = ['top15',
                    'classification_bar',
                    'classification_bar_w_uc',
                    'ecdf',
                    'tpr_fps_curve',
                    'pr_curve',
                    'l1_heatmap'
                    ]
    table_names = ['metric_table',
                   'number_table'
                   ]
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
            template_vars[table + '_' + organism] = pd.read_table(table_path, index_col=0, sep='\t')\
                .to_html(float_format='%.5g')

    html_out = template.render(template_vars)
    with open(output_path / f'generated_report_{title}.html', 'w+') as f:
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
        parser.add_argument('--input',
                            help="Output from SnakeMake flow.\n"
                                 "Format is [name]/ouput/reads.",
                            required=True
                            )
        parser.add_argument('--output',
                            help='Location of report/figures',
                            required=True
                            )
        parser.add_argument('--classifier',
                            choices={'bracken', 'desamba', 'mash', 'minimap2', 'kraken2', 'kma',
                                     'metaphlan', 'centrifuge', 'mmseqs2', 'motus'},
                            nargs='+',
                            default='all',
                            help='Which classifier(s) to use [all]',
                            metavar='CLASSIFIER')

        return parser.parse_args(args)

    def __get_config_data(self) -> Dict[str, Any]:
        """
        Get the snakemake configuration data.
        :return: Configuration data
        """
        return {
            'output': self._args.output,
            'input': self._args.input,
            'classifier': self._args.classifier,
        }

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        generate_plots(Path(self._args.output), Path(self._args.input), self._args.classifier, self.__get_config_data())
        logger.info(f"Images created")

        title = Path(self._args.input).parents[0].stem

        generate_report(Path(self._args.output), title)


if __name__ == '__main__':
    logger = logging.getLogger('amr_report')
    logger.setLevel(logging.INFO)

    logger.addHandler(RichHandler(show_path=False, level=logging.INFO))

    report = Report()
    report.run()
