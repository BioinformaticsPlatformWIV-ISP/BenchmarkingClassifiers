import argparse
import logging
import yaml
from pathlib import Path
from rich.logging import RichHandler
from typing import Dict, Any, Optional, Sequence
from datetime import datetime

from utils.command import Command


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

    # Create string with SnakeMake params
    snakemake_params = [
        f'snakemake',
        f"--snakefile {snakefile}",
        f'--configfile {config_path}',
        f'--cores',
        f'-p',
        f'--stats {dir_working / ("stats_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + ".json")}',
        f'-r',
        f'-n'
    ]
    # Execute  snakemake
    if dry_run:
        command_ = Command(' '.join(snakemake_params))
    else:
        command_ = Command(' '.join(snakemake_params[:-1]))

    command_.run_command(dir_working)

    # Check command output
    if command_.exit_code != 0:
        raise RuntimeError(f"Error executing snakemake workflow:\n{command_.stderr}")
    with (dir_working / 'stderr.log').open('w') as handle:
        handle.write(command_.stderr)


class Classifying(object):
    """
    This class contains is used to perform the classifying .
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the an classifying context instance.
        """
        self._args = Classifying._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('--input',
                            help="FASTA/Q input file",
                            required=True,
                            type=lambda x: Path(x).resolve()
                            )
        parser.add_argument('--output',
                            help='Output of snakemake',
                            required=True,
                            type=lambda x: Path(x).resolve())
        parser.add_argument('--dir_working',
                            help='Working directory. If not set, it will be the output path.',
                            default=None
                            )
        parser.add_argument('--config_file',
                            help='YAML file with full path to DBs',
                            default=(Path(__file__).parent / 'config' / 'db.yml').resolve()
                            )

        parser.add_argument('--dry_run', help="Run the SnakeMake file in Dry Run mode", action='store_true')
        parser.add_argument('--no_filter', help='Disable filtering of reads', action='store_true')
        parser.add_argument('--classifier',
                            choices={'bracken', 'ccmetagen', 'centrifuge', 'kaiju', 'kma', 'kraken2', 'mash', 'metaphlan',
                                     'minimap2', 'mmseqs2', 'motus'},
                            nargs='+',
                            default={'bracken', 'ccmetagen', 'centrifuge', 'kaiju', 'kma', 'kraken2', 'metaphlan',
                                     'mmseqs2', 'motus'},
                            help='Which classifier(s) to use. Choices are %(choices)s; space separated.\n '
                                 'Default are all.',
                            metavar=''
                            )
        parser.add_argument('--truth',
                            help='Location to the ground truth file. This will update the taxonomy in the file to the taxdump used in the '
                                 'classifiers.',
                            required=True
                            )

        parser_args = parser.parse_args(args)

        # Set working directory to output directory if not specified
        if parser_args.dir_working is None:
            parser_args.dir_working = parser_args.output / 'snakemake_conf'

        return parser_args

    def _check_input(self):
        """
        Check if input file really exists
        """
        if not self._args.input.exists():
            logger.error(f'Input file {self._args.input} does not exist!')
            exit()

    def __get_config_data(self) -> Dict[str, Any]:
        """
        Get the snakemake configuration data.
        :return: Configuration data
        """

        with open(self._args.config_file, 'r') as f:
            dbs = yaml.load(f, Loader=yaml.FullLoader)

        return {
            'output': str(self._args.output),
            'input': {'fastq': str(self._args.input)},
            'no_filter': self._args.no_filter,
            'db': dbs,
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
        self._check_input()

        run_snakemake(
            self._args.dir_working,
            (Path(__file__).parent / 'smk' / 'map_classifiers.smk').resolve(),
            self._args.dry_run,
            self.__get_config_data())


if __name__ == '__main__':

    logger = logging.getLogger()
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
