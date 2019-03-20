#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
#from matplotlib import pyplot, mpl, gridspec
from matplotlib import pyplot, gridspec
import matplotlib as mpl
from pylab import *
from scipy import *
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Bbox

parser = OptionParser(usage="usage: %prog [options] <taxonomy_heatmap.txt>")

parser.add_option('--output', '-o', action="store", default="blast_map", type='string', dest="output",
                  help='Parsed blast output file header (default=%default).  Will write (output).png.')
parser.add_option('--max', '-m', action="store", default=1000, type='int', dest="max",
                  help='Maximum number of sequences to plot (default=%default). The longest contigs up to this number will be plotted.')
parser.add_option('--interesting', '-i', action="store_true", default=None, dest="interesting",
                  help='Only plot sequences which have more than one different taxonomic hit. (default=%default).')
parser.add_option('--grouping', '-g', action="store", type='string', default=None, dest="group_file",
                  help='A file describing how to group contigs.  If a file with the extention .agp is given, the file will be processed like an AGP.  Otherwise, a two column tab-delimited file will be expected with the unit in the first column and the group id in the second column')
parser.add_option('--no_scaling', action="store_true", default=False, dest="no_scaling",
                  help='Depreciated.  Scaling is no longer an option.')

(options, args) = parser.parse_args()

if len(args) < 1:
    parser.error("Must supply taxonomy heatmap and agp file output.")
    sys.exit(-1)

class HeatmapUnit:
    """This class represents a single sequence unit's taxonomic deciles"""
    def __init__(self,tax_string):
        self.id = tax_string.pop(0)
        self.length = tax_string.pop(0)
        self.decile_orgs = []
        self.decile_alphas = []
        self.orgs = []
        #self.colors = []

        for decile in tax_string:
            orgnum,alpha = decile.split(".")
            self.decile_orgs.append(orgnum)
            self.decile_alphas.append(alpha)
            #self.colors.append(mcolors.to_rgba(orgnum,alpha=alpha))
            if ((not orgnum in self.orgs)):
                #print orgnum, "not seen yet.  Incrementing"
                self.orgs.append(orgnum)

    def num_orgs(self):
        return len(self.orgs)

def __parse_input():
    orgnames = {}
    data = {}
    
    for line in open(args[0], 'r'):
        if (line.startswith("#")):
            key, values = line.rstrip().split("  ")
            keys = []
            keys = values.split(";")
            for set in keys:
                key,value = set.split("=")
                value = int(value)
                orgnames[value]=key
        else:
            hmu = HeatmapUnit(line.rstrip().split("\t"))
            if (not options.interesting or hmu.num_orgs() > 1):
                data[hmu.id] = hmu
                #print "Plotting", hmu.id
                #print hmu.orgs
                #print hmu.decile_orgs
                #print hmu.decile_alphas
            else:
                print "Skipping", hmu.id, "since it only has", hmu.num_orgs(), "unique taxanomic hits."
                
    return orgnames,data

def __plot_data(orgnames,data,output,max_length,groups):

    cmap = cm.spectral
    mcolors = cm.ScalarMappable()
    mcolors.set_cmap(cmap)
    marray = np.array(orgnames.keys())
    mcolors.set_array(np.array(orgnames.keys()))
    mcolors.autoscale()

    nunits = min([len(data),options.max])    
    print "Creating figure for", nunits, "plotted units..."

    #Creating figure
    fig_width = 10
    fig_length = (nunits/5)+3
    if fig_length < 4:
        print "Figure height too small, increasing to 4 inches."
        fig_length = 4;
    fig = figure(figsize=(fig_width,fig_length),facecolor=None,frameon=True)
    fig.text(0.5,1-0.5/fig_length,"Taxonomic Hits", ha='center', va='top',fontsize=24)
    #subplots_adjust(left = 1/fig_width, bottom = 3/fig_length, right = 1 - 1/fig_width, top = 1 - 1/fig_length,wspace=0,hspace=0.1)
    subplots_adjust(left = 0.5/fig_width, bottom = 2.0/fig_length, right = 1 - 0.5/fig_width, top = 1 - 1.0/fig_length,wspace=0,hspace=0.1)
    print "Created",fig_width,"x",fig_length,"inch image."

    plottedorgs = []
    print "Plotting units..."
    group_row = 0
    unit_row = 0

    for group in (sorted(groups.keys())):
        group_size = 0

        if (unit_row+1 > options.max):
            print "WANRING! Stopped after", options.max,"units.  To print more units please use the -m option, or control the specific units you want to plot with -g."
            break
        for unit in (groups[group]):
            group_size += 1
            tax_data = data[unit]
            #print tax_data
            colors = []
            for org,alpha in zip(tax_data.decile_orgs, tax_data.decile_alphas):
                if alpha == "00":
                    alpha = 100
                if org == "0":
                    alpha = 0
                alpha = float(alpha)/100.0
                colors.append(mcolors.to_rgba(org,alpha=alpha))
                #print org, mcolors.to_rgba(org,alpha=alpha)

            if (unit_row+1 > options.max):
                continue

            #print "Plotting",unit
            mydata = []
            mydata.append(colors)
            
            #plot for bar
            ax = subplot2grid((nunits,7),(unit_row,1),colspan=5)
            ax.tick_params(axis='both',bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')
            imshow(mydata, aspect='auto', origin='lower')

            #Adjust Bar size if scaled
            pos = list(ax.get_position().bounds)

            #Print Text label for unit
            ax.text(0.5, 0.0, tax_data.id, fontsize=8, horizontalalignment='left', va='center')

            #Print Size label for unit
            kb = float(tax_data.length)/1000

            #Add printed tax data to 
            for x in tax_data.orgs:
                if not x in plottedorgs:
                    plottedorgs.append(x)

            #plot for size bar
            if (unit_row > 0):
                size_ax = subplot2grid((nunits,7),(unit_row,6),sharex=first_ax)
            else:
                size_ax = subplot2grid((nunits,7),(unit_row,6))
                first_ax = size_ax
            size_ax.tick_params(axis='both',bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')
            size_ax.text(0.5,0.5,"%.1f"%round(kb)+"kb", fontsize=8, color="black", ha='left', va='center')
            size_ax.set_axis_off()
            barh(0,kb,height=1.0,color='0.75',linewidth=0)

            #increment for next unit
            unit_row += 1

            
        #plot for group
        group_ax = subplot2grid((nunits,7),(group_row,0),rowspan=group_size)
        group_ax.tick_params(axis='both',bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')
        group_ax.text(0.01,1.0, group, fontsize=8, ha='left', va='top')
        group_row = group_row + group_size

    print "Printed", min(nunits,options.max), "plots."

    llabels = []
    lpoints = []

    for value in plottedorgs:
        #print value, orgnames[int(value)], mcolors.to_rgba(value)
        color = mcolors.to_rgba(value)
        if value == "0":
            continue
        p = Rectangle((0, 0), 1, 1, fc=color)
        lpoints.append(p)
        llabels.append(orgnames[int(value)])

    fp = FontProperties(size='small')
    fig.legend(lpoints,llabels,bbox_to_anchor=(-0.2,0,1.4,.102),frameon=1,loc=8,labelspacing=0.01,ncol=3,markerscale=0.5,prop=fp)
    figure_name = output + ".blast_map.png"
    print "Creating", figure_name
    savefig(figure_name)
    return 0

def main():

    orgnames, data, = __parse_input()
    groups = {}
    max_length = 0
    units = []
    
    for id,hmu in sorted(data.iteritems(), key=lambda item: int(item[1].length),reverse=True ):
        #print id,hmu.length,max_length
        units.append(id)
        if int(hmu.length) > int(max_length):
            max_length = hmu.length

    print "Total Units:", len(data) 
    print "Longest unit is", max_length

    if (options.group_file):
        print "Grouping using", options.group_file
        agp_re = re.compile(".*\.agp")
        match = agp_re.match(options.group_file)
        if (match):
            for line in open(options.group_file, 'r'):
                agp_line = line.rstrip().split("\t")
                if agp_line[0] not in groups.keys():
                    groups[agp_line[0]] = []
                if (agp_line[4] == "W"):
                    groups[agp_line[0]].append(agp_line[5])
        else:
            print "Need to implement non-agp group file."
            sys.exit(-1)
    else:
        groups["no groups"] = units

    return __plot_data(orgnames,data,options.output,max_length,groups)

if __name__ == "__main__":
    sys.exit(main())
