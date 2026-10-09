"""Allow ``python -m pufprobe`` as an alias for the console script."""

from pufprobe.cli import cli

if __name__ == "__main__":
    cli()
