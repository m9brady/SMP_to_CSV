import os

import click
import matplotlib.pyplot as plt
from snowmicrotoolset import SMP
from snowmicrotoolset import __version__ as pkg_version
from tqdm import tqdm


@click.command(
    context_settings={
        "help_option_names": ["-h", "--help"],
        "max_content_width": 120,
    },
)
@click.option(
    "--input_folder",
    "-i",
    #default="./indata",
    #show_default=True,
    help="Location of input .smp files"
)
@click.option(
    "--output_folder",
    "-o",
    default="./outdata",
    show_default=True,
    help="Location of output csv/png files"
)
@click.pass_context
@click.version_option(version=pkg_version, prog_name="snowmicrotoolset")
def cli(ctx, input_folder, output_folder):
    '''
    Convert SLF .PNT files to CSV/PNG formats
    '''
    if input_folder is None:
        click.echo(
            "\n***************************\n" +\
            "** No input folder given **\n" +\
            "***************************\n"
        )
        click.echo(ctx.get_help())
        ctx.exit()
    if not os.path.exists(input_folder):
        raise IOError("No such folder: %s" % input_folder)
    if not os.path.exists(output_folder):
        try:
            os.makedirs(output_folder)
        except (IOError, OSError):
            raise IOError("Cannot create output folder: %s" % output_folder)
    # walk through subdirectories of input_data, looking for SMP .pnt files
    to_process = [
        os.path.join(root, f) 
        for root, folders, files in os.walk(input_folder)
        for f in files
        if f.endswith(".pnt")
    ]
    if len(to_process) == 0:
        raise IOError("No PNT files found in folder: %s" % input_folder)
    # process each file sequentially
    for pnt in tqdm(to_process):
        file_name = os.path.basename(pnt)
        unique_id = os.path.splitext(file_name)[0]
        p = SMP(pnt)
        date = "%s%s%s" % (
            p.header["Year"],
            str(p.header["Month"]).zfill(2),
            str(p.header["Day"]).zfill(2)
        )
        csv_name = os.path.join(
            output_folder,
            "SMP_%s_%s.csv" % (date, unique_id)
        )
        png_name = csv_name.rstrip('.csv') + '.png'
        if not os.path.isfile(csv_name):
            p.export_to_csv(csv_name)
        # ensure matplotlib interactive is off
        if plt.isinteractive():
            plt.ioff()
        if not os.path.isfile(png_name):
            p.plot_quicklook(png_name)
