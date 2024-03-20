import argparse
import logging
import sys

import yaml
import snakemake
from pathlib import Path
from rich.logging import RichHandler
from typing import Dict, Any, Optional, Sequence
from datetime import datetime


def run_snakemake(dir_working: Path, snakefile: Path, dry_run: bool, config_data: Dict[str, Any]) -> None:
    """
    Runs snakemake.
    :param dir_working: Directory in which the *.smk will be executed
    :param snakefile: Snakefile
    :param dry_run: If True, will execute the SnakeMake command with parameter -n
    :param config_data: Snakemake configuration data
    :return:None
    """
    # Create working directory
    if not dir_working.exists():
        dir_working.mkdir(parents=True)

    # Create snakemake config
    config_path = (dir_working / 'snakemake_config.yml').resolve()
    with config_path.open('w') as handle:
        yaml.dump(config_data, handle)
    logger.info(f"Snakemake config created: {config_path}")

    snakemake.snakemake(snakefile=snakefile,
                        configfiles=[config_path],
                        workdir=dir_working,
                        stats=dir_working / ("stats_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + ".json"),
                        printreason=True,
                        printshellcmds=True,
                        cores=32,
                        dryrun=dry_run
                        )


class Classifying(object):
    """
    This class contains is used to perform the classifying .
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the classifying context instance.
        """
        self._args = Classifying._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('--input',
                            help="FASTA/Q input file",
                            required=True,
                            type=lambda x: Path(x).resolve()
                            )
        parser.add_argument('--truth',
                            help='Location to the ground truth file.',
                            required=True
                            )

        parser.add_argument('--output',
                            help='Output of snakemake',
                            required=True,
                            type=lambda x: Path(x).resolve())
        parser.add_argument('--config_file',
                            help='YAML file with location of tools and database.\n'
                                 'Template can be found under "/snakemake/config/classifiers.yml.template"',
                            required=True
                            )
        parser.add_argument('--classifier',
                            nargs='+',
                            default=['bracken', 'ccmetagen', 'centrifuge', 'kaiju', 'kma', 'kraken2', 'metaphlan',
                                     'mmseqs2', 'motus'],
                            help='Which classifier(s) to use. \n'
                                 'Default choices are %(default)s; space separated.\n'
                                 'If not set, all default choices will be used.',
                            metavar=''
                            )
        parser.add_argument('--dir_working',
                            help='Working directory. \n'
                                 'If not set, it will be {output path} + "snakemake_conf".',
                            default=None
                            )
        parser.add_argument('--dry_run', help="Run the SnakeMake file in Dry Run mode", action='store_true')
        parser.add_argument('--no_filter', help='Disable filtering of reads (length > 1000 and quality > 7)', action='store_true')

        parser_args = parser.parse_args(args)

        return parser_args

    def _check_input_arguments(self):
        """
        Check if input file really exists and arguments make sense
        """
        if not self._args.input.exists():
            logger.error(f'Input file {self._args.input} does not exist!')
            sys.exit(1)

        mutual_inclusive = {'ccmetagen': 'kma', 'bracken': 'kraken2'}
        for name, dependency in mutual_inclusive.items():
            if name in self._args.classifier and dependency not in self._args.classifier:
                logger.error(f"Error: '{dependency}' is required when '{name}' is selected.")
                sys.exit(1)

        # Set working directory to output directory if not specified
        if self._args.dir_working is None:
            self._args.dir_working = self._args.output / 'snakemake_conf'
            logger.info(f"Working DIR set to {self._args.dir_working}")

    def __get_config_data(self) -> Dict[str, Any]:
        """
        Get the snakemake configuration data.
        :return: Configuration data
        """

        with open(self._args.config_file, 'r') as f:
            locations = yaml.load(f, Loader=yaml.FullLoader)

        return {
            'output': str(self._args.output),
            'input': {'fastq': str(self._args.input)},
            'no_filter': self._args.no_filter,
            'locations': locations,
            'classifiers': self._args.classifier,
            'truth_file': self._args.truth,
            'python_scripts': str(Path(__file__).parent.resolve() / 'scripts')
        }

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Check if input exists
        self._check_input_arguments()

        run_snakemake(
            self._args.dir_working,
            (Path(__file__).parent / 'smk' / 'map_classifiers.smk').resolve(),
            self._args.dry_run,
            self.__get_config_data())


if __name__ == '__main__':

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    if not logger.handlers:
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        # create formatter
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        # add formatter to ch
        ch.setFormatter(formatter)
        # add ch to logger
        # logger.addHandler(ch)
        logger.addHandler(RichHandler(show_path=False, level=1))

    classify = Classifying()
    classify.run()
