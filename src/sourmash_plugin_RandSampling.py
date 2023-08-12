"""\
    RandSampling plugin performs wighted-random-sampling to extract a max number of k-mers from a signature using their abundance distribution.
    
    The plugin main aim is to extract kmers from large signatures with 
    different coverages (abundances) to create a smaller signature with 
    a more uniform coverage distribution.
    
    The --force option allows to force the extraction of a signature in case its empty or has less than max_kmers.
    This is useful in workflows (e.g. snakemake). 
"""

usage="""
   sourmash scripts randsampling
"""

epilog="""
See https://github.com/sourmash_plugin_RandSampling for more examples.

Need help? Have questions? Ask at https://github.com/sourmash_plugin_RandSampling/issues OR http://github.com/sourmash/issues!
"""

import argparse
import sourmash

from sourmash.index import LinearIndex
from sourmash.logging import debug_literal
from sourmash.plugins import CommandLinePlugin
from sourmash.save_load import (Base_SaveSignaturesToLocation, _get_signatures_from_rust)


from sourmash import sourmash_args
from sourmash.logging import debug_literal, error, notify
from sourmash.plugins import CommandLinePlugin
from sourmash.cli.utils import (add_ksize_arg, add_moltype_args)

import shutil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys


def plot_abundances(sig, selected_abundances, output_sig, show = False):
    def list_to_hist(l):
        hist = {}
        for i in l:
            if i in hist:
                hist[i] += 1
            else:
                hist[i] = 1
        return hist

    custom_median = lambda numbers: sorted(numbers)[len(numbers)//2 - 1] if len(numbers) % 2 == 0 else sorted(numbers)[len(numbers)//2]


    sig_minhash = sig.minhash
    sig_hashes = sig_minhash.hashes
    abund_to_count = list_to_hist(sig_hashes.values())
    abund_to_count_df = pd.DataFrame(list(abund_to_count.items()), columns=['Abundance', 'Count'])
    abund_to_count_df = abund_to_count_df.sort_values('Abundance')
    median_count = custom_median(abund_to_count_df['Count']) # np.median(abund_to_count_df['Count'])

    # Set Seaborn style
    sns.set_style("whitegrid")
    sns.set_context("notebook", font_scale=1.2)

    plt.figure(figsize=(10, 6))

    sns.lineplot(x='Abundance', y='Count', data=abund_to_count_df, lw=2, color='blue')
    # sns.lineplot(x=abund_to_count_df['Abundance'], y=abund_to_count_df['Count'], data=abund_to_count_df, lw=2, color='blue')

    plt.axhline(y=median_count, color='r', linestyle='--', lw=2, label="Median Abundance Frequency")
    print(f"Median Count: {median_count}")
    # get abundance at median count
    abund_to_count_df.to_csv("abund_to_count_df.csv")
    
    median_abundance = abund_to_count_df['Abundance'][abund_to_count_df['Count'] == median_count].iloc[0]
    print(f"Median Abundance: {median_abundance}")
    
    plt.scatter([median_abundance], [median_count], color='red', s=100)
    
    selected_counts = [abund_to_count[abundance] for abundance in selected_abundances]
    plt.scatter(selected_abundances, selected_counts, color='magenta', edgecolor='black', linewidth=1, s=10, label="Selected Abundances (k-mers)")
    plt.xlabel('Abundance')
    plt.ylabel('Frequence')
    plt.title('Abundance vs Frequency Distribution')
    plt.legend()

    plt.yscale('log')
    plt.xscale('log')

    # Saving the plot as a high-resolution image
    
    # get extensoin from output_sig
    # if output_sig has extension
    if "." in output_sig:
        _extension = output_sig.split(".")[-1]
        output_image = output_sig.replace(_extension, "png")
    else:
        output_image = output_sig + ".png"

    notify(f"Saving plot to {output_image}")
    plt.savefig(output_image, dpi=500)
    
    if show:
        plt.show()



class Command_RandSampling(CommandLinePlugin):
    command = 'randsampling'             # 'scripts <command>'
    description = __doc__       # output with -h
    usage = usage               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline



    def __init__(self, subparser):
        super().__init__(subparser)
        
        subparser.add_argument('--sig', type=str, help='path to signature file', required=True)
        subparser.add_argument('-m', '--max-kmers', type=int, help='maximum number of kmers to extract', required=True)
        subparser.add_argument('-k', type=int, help='kmer size', required=True)
        subparser.add_argument('-o', '--out', type=str, help='path to output signature file', required=True)
        subparser.add_argument('--plot', action='store_true', help='plot abundance distribution', required=False, default=False)
        subparser.add_argument('--force', action='store_true', help='force write new signature if empty or <max_kmers', required=False, default=True)

        debug_literal('RUNNING cmd randsampling __init__')       



    def main(self, args):
        super().main(args)

        sig = sourmash.load_one_signature(args.sig, ksize=args.k)
        minhash = sig.minhash


        if len(minhash.hashes) == 0:
            notify("Empty signature detected")
            if args.force:
                notify(f"Force copying {args.sig} to {args.out}")
                with sourmash.sourmash_args.FileOutput(args.out, 'wt') as fp:
                    sourmash.save_signatures([sig], fp=fp)
                sys.exit(0)
            else:
                notify("Exiting.")
                sys.exit(1)

        # Check if there's abundance in the minhash
        if not minhash.track_abundance:
            error("There is no abundance in the minhash.")
            sys.exit(1)


        if len(minhash.hashes) <= args.max_kmers:
            print(f"Number of kmers in the signature is less than or equal to the predefined maximum number of kmers ({args.max_kmers}).")
            if args.force:
                print(f"Force copying {args.sig} to {args.out}")
                shutil.copyfile(args.sig, args.out)
                sys.exit(0)
            else:
                print("Exiting.")
                sys.exit(1)                


        ##################################
        # Perform weighted subsampling
        ##################################

        kmers = list(minhash.hashes.keys())
        abundance = list(minhash.hashes.values())

        total_abundance = sum(abundance)

        # Normalize the abundance values to create a probability distribution
        probabilities = [a / total_abundance for a in abundance]

        # If due to numerical precision, probabilities do not sum to 1, fix by scaling them
        if sum(probabilities) != 1.0:
            scaling_factor = 1.0 / sum(probabilities)
            probabilities = [p * scaling_factor for p in probabilities]

        # Perform weighted random subsampling without replacement up to max_kmers
        selected_kmers = np.random.choice(kmers, size=min(args.max_kmers, len(kmers)), p=probabilities, replace=False)

        ##################################

        notify(f"Extracted {len(selected_kmers)} kmers out of {len(minhash.hashes)} kmers.")

        final_mh = sig.minhash.copy_and_clear().flatten()
        final_mh.add_many(selected_kmers)
        final_mh = final_mh.inflate(sig.minhash)


        finalSig = sourmash.SourmashSignature(final_mh)
        with sourmash.sourmash_args.FileOutput(args.out, 'wt') as fp:
            sourmash.save_signatures([finalSig], fp=fp)

        if args.plot:
            final_mh_abundances = final_mh.hashes.values()
            plot_abundances(sig, final_mh_abundances, args.out, False)
        else:
            notify("Plotting is disabled. Rerun with --plot to generate a plot of abundance distribution.")
            
        
        notify(f"Saved the subsampled signature to {args.out}")
        