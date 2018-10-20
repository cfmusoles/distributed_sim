#!/usr/bin/env python
# Create plots from Apprentice2 mosaic trace file or pat_report ap2-xml file
# 
# Copyright Cray UK Ltd, 
# Harvey Richardson and Michael Bareford, Copyright 2017
# v1.0


import sys
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons
from matplotlib.widgets import Slider
import matplotlib.patches as mpatches


def have_package_pylab():
  try:
    import pylab
    return True
  except ImportError:
    return False

  
def help():
  print """usage: app2mosaic_full_display file

  Arguments:

    -d id [-point-size ps]     Choose dataset (and plot point size)
    -node [-inter-node-only]   Display node-to-node mosaic
    -node-summary              Display summary of off/on-node communications
    -display-output file       Filename used to output plot image
    file                       ap2 xml input file

 Description

  This script parses output from perftools to obtain the P2P communication data.
  The input file is generated via pat_report.

    pat_report -f ap2-xml experiment.ap2

  You could then run this script on the output.

    app2mosaic_full_display -d 0 experiment.ap2-xml

  A rank-to-rank mosaic is presented for the dataset indicated by "-d 0".
  A node-to-node display is shown if the "-node" argument is present.
  The "-node" argument can be accompanied by the "-inter-node-only" qualifier,
  which means the plot will only show inter node comms data and exclude any
  communications between PEs on the same node.
  The "-node-summary" argument displays a stacked bar plot for all the different
  datasets that shows which proportion of the communications was on/off-node.

"""

def parse_args():
  global ifn
  global idata
  global node2node
  global nodesummary
  global internodeonly
  global display_output
  global user_pt_size

  for i in range(1,len(sys.argv)):
    arg=sys.argv[i]
    if arg == "-d":
      i += 1
      idata = int(sys.argv[i])
    elif arg == "-point-size":
      i += 1
      user_pt_size = float(sys.argv[i])
    elif arg == "-node":
      node2node = True
      i += 1
      internodeonly = ("-inter-node-only" == sys.argv[i].lower())
    elif arg == "-node-summary":
      node2node = True
      nodesummary = True
    elif arg == "-display-output":
      i += 1
      display_output = sys.argv[i]
    elif arg == "-h" or arg == "-help":
      help()
      sys.exit();
    else:
      ifn = arg
    
  return

  
    
def parse_pat_data():
  global ifn
  global nranks, nnodes
  global pe_node_xyz
  global pe_node_id
  global node_ids
  global data_dict
  
  re_pes = re.compile(" pes='(\d+)'")
  re_data_table = re.compile("<data_table name='(\S+P2P\S+)'")
  re_coll_data_table = re.compile("<data_table name='(\S+Coll\S+)'")
  re_data_table_end = re.compile("</data_table")
  re_data_tag = re.compile(".*data")
  re_data = re.compile("<data ctidx='(\d+)' first_pe='(\d+)'")
  re_data_end = re.compile("</data")
  re_pdata_node_id = re.compile("<pdata name='node_id' first_pe='(\d+)'")
  re_pdata_mesh_xyz = re.compile("<pdata name='mesh_xyz' first_pe='(\d+)'")
  re_pdata_end = re.compile("</pdata")

  inside_P2P_table = False
  data_label = ""

  with open(ifn) as fin:
    
    for line in fin:

      if "cpat_version" in line:
        pat_version_string = line[line.find("'")+1:-2]
        print pat_version_string
      
      m = re_pes.match(line)
      if m:
        nranks = int(m.group(1))
        data_2d = np.zeros((nranks,nranks))
        continue

      m = re_data_table.match(line)
      if m:
        inside_P2P_table = True
        data_label = m.group(1)

      """
      m = re_coll_data_table.match(line)
      if m:
        inside_P2P_table = True
        data_label = m.group(1)
      """
        
      # Get mesh_xyz
      # pdata name='mesh_xyz' first_pe='<pe>'
      # Aries (X optical, Y group, Z backplane)
      #  num values (or *n values)
      #  values are of form 0xXXXXYYYYZZZZ no leading zeros
      #   So  0x200030002  is 2,3,2
      m = re_pdata_mesh_xyz.match(line)
      if m:
        first_pe = int(m.group(1))
        j = first_pe
        row_data = []
        for line in fin:
          if re_pdata_end.match(line):
            break
          else:
            row_data += line.rstrip().split(' ')
        repeat = 1
        for val in row_data[1::]:
          if val[0] == '*':
            repeat = int(val[1:])
            continue
          v = int(val,16)
          for ij in range(repeat):
            pe_node_xyz[j] = v
            j += 1
          repeat = 1

      # Get pdata node_id info for node mapping
      # Put if nodal around following
      m = re_pdata_node_id.match(line)
      if m:
        first_pe = int(m.group(1))
        j = first_pe
        # line+next rows have the actual data 
        # npes then data per pe with possible *n v repeat counts
        # note the *n v can be split over line boundaries
        row_data = []
        for line in fin:
          if re_pdata_end.match(line):
            break
          else:
            row_data += line.rstrip().split(' ')
        repeat = 1
        for val in row_data[1::]:
          if val[0] == '*':
            repeat = int(val[1:])
            continue
          v = int(val)
          for ij in range(repeat):
            pe_node_id[j] = v
            if v not in node_ids:
              node_ids.append(v)
            j += 1
          repeat = 1
        
      if re_data_table_end.match(line):
        if len(data_label) > 0:
          data_dict[data_label] = data_2d
          data_2d = np.zeros((nranks,nranks))
          data_label = ""
        inside_P2P_table = False

      if inside_P2P_table:
        if re_data_tag.match(line):
          m = re_data.match(line)
          if m:
            pe_idx = int(m.group(1))
            j = int(m.group(2))
        else:
          # line+next rows have the actual data 
          # npes then data per pe with possible *n v repeat counts
          row_data = line.rstrip().split(' ')
          for line in fin:
            if re_data_end.match(line):
              break
            else:
              row_data += line.rstrip().split(' ')
              
          repeat = 1
          for val in row_data[1::]:
            if val[0] == '*':
              repeat = int(val[1:])
              continue
            v = int(val)
            v = float(val)
            for ij in range(repeat):
              data_2d[pe_idx,j] = v
              j += 1
            repeat = 1

  nnodes = len(node_ids)
  return


def dist_xyz(a, b):
  dx = abs ( a/0x100000000 - b/0x100000000 )
  dy = abs ( a/0x10000%0x10000 -b/0x10000%0x10000 )
  dz = abs ( a%0x10000 - b%0x10000 )
  return dx, dy, dz

def calculate_network_distances():
  global nranks
  global pe_node_xyz
  global pe_node_id
  global data_dict
  
  data_2d = np.zeros((nranks,nranks))

  counts_2d = data_dict["Data_P2P_Count"]
  
  for i in range(nranks):
    for j in range(nranks):
      if counts_2d[i,j] == 0:
        continue
      
      dx,dy,dz = dist_xyz(pe_node_xyz[i], pe_node_xyz[j])
      if dx > 0:
        # pes communicate over rank 3 network (optical)
        data_2d[i,j] = 4
      elif dy > 0:
        # pes communicate over rank 2 network (copper/group)
        data_2d[i,j] = 3
      elif dz > 0:
        # pes communicate over rank 1 network (backplane)
        data_2d[i,j] = 2
      elif pe_node_id[i] != pe_node_id[j]:
        # pes are on same blade
        data_2d[i,j] = 2
      else:
        # pes are on same node
        data_2d[i,j] = 1
        
  data_dict["Data_P2P_NetworkDist"] = data_2d


def collate_node_data():
  global nranks, nnodes
  global pe_node_xyz
  global pe_node_id
  global data_dict
  global node_data_dict
  global node_summary_dict

  node_data_dict = {}
  
  for key in data_dict:
    n2n_data = np.zeros((nnodes,nnodes))
    p2p_data = data_dict[key]
    min_key = "min" in key.lower()
    max_key = "max" in key.lower()
    network_key = "network" in key.lower()
    
    for i in range(nranks):
      for j in range(nranks):
        if 0 == p2p_data[i,j]:
          continue
        
        inode = pe_node_id[i]
        jnode = pe_node_id[j]

        if inode==jnode and internodeonly:
          continue
        
        ind = node_ids.index(inode)
        jnd = node_ids.index(jnode)

        if min_key:
          if 0 == n2n_data[ind, jnd]:
            n2n_data[ind, jnd] = p2p_data[i,j]
          elif p2p_data[i,j] < n2n_data[ind, jnd]:
            n2n_data[ind, jnd] = p2p_data[i,j]
        elif max_key:
          if p2p_data[i,j] > n2n_data[ind, jnd]:
            n2n_data[ind, jnd] = p2p_data[i,j]
        elif network_key:
          n2n_data[ind, jnd] = p2p_data[i,j]
        else:    
          n2n_data[ind, jnd] += p2p_data[i,j]
        
    node_data_dict[key] = n2n_data
          
          
        
def get_simple_data_label(label):
  global data_label_prefix
  return label[len(data_label_prefix):]

def get_snappy_data_label(label):
  global data_label_prefix
  lb = label[len(data_label_prefix):]

  if lb == "SendTime":
    lb = "Send Time"
  elif lb == "SendTimeMax":
    lb = "Maximum Send Time"
  elif lb == "SendTimeMin":
    lb = "Minimum Send Time"
  elif lb == "RecvTime":
    lb = "Receive Time"
  elif lb == "RecvTimeMax":
    lb = "Maximum Receive Time"
  elif lb == "RecvTimeMin":
    lb = "Minimum Receive Time"
  elif lb == "NetworkDist":
    lb = "Network Separation"

  return lb
    

def prepare_plot_data(label):
  global nranks, nnnodes
  global node2node
  global data_dict, node_data_dict
  global data_label
  global data_2d
  global data_unit_label
  global data_min, data_max
  global data_all

  data_label = data_label_prefix + label
  if node2node:
    data_2d = node_data_dict[data_label]
  else:
    data_2d = data_dict[data_label]

  data_unit_label = ""
  if "time" in data_label.lower():
    # convert times from nanoseconds to seconds
    data_2d /= 1.0e9
    data_unit_label = "s"
  elif "bytes" in data_label.lower():
    # convert bytes to megabytes
    data_2d /= 1.0e6
    data_unit_label = "MB"

  data_1d = data_2d.flatten()
  data_min = np.min(data_1d)
  data_max = np.max(data_1d)
  
  data_all = {}
  data_all["cnt"] = len(data_1d[data_1d>0.0])
  data_all["xpts"] = np.zeros(data_all["cnt"])
  data_all["ypts"] = np.zeros(data_all["cnt"])
  data_all["vals"] = np.zeros(data_all["cnt"])  

  nitems = nnodes if node2node else nranks
  k = 0
  for i in range(nitems):
    for j in range(nitems):
      if data_2d[i,j] > 0.0:
        data_all["xpts"][k] = i
        data_all["ypts"][k] = j
        data_all["vals"][k] = data_2d[i,j]
        k += 1
    

def prepare_summary_data():
  global internodeonly
  global nranks
  global data_dict
  global network_labels
  global network_summary_dict
  
  if "Data_P2P_NetworkDist" not in data_dict:
    print "Error, network info missing."
    return

  network_data = data_dict["Data_P2P_NetworkDist"]
  network_summary_dict = {}
  
  for key in data_dict:
    if "network" in key.lower():
      continue

    if "min" in key.lower() or "max" in key.lower():
      continue

    summary_data = np.zeros(len(network_labels))
    p2p_data = data_dict[key]
    for i in range(nranks):
      for j in range(nranks):
        summary_data[int(network_data[i,j])-1] += p2p_data[i,j]

    # normalise summary totals
    summary_data_sum = np.sum(summary_data)
    summary_data /= summary_data_sum
    network_summary_dict[key] = summary_data

  # calculate baseline summary, i.e., the proportion
  # of communicating ranks that are on-node, backplane and so on.
  summary_data = np.zeros(len(network_labels))
  count_data = data_dict["Data_P2P_Count"]
  for i in range(nranks):
    for j in range(nranks):
      if 0 < count_data[i,j]:
        summary_data[int(network_data[i,j])-1] += 1
  summary_data_sum = np.sum(summary_data)
  summary_data /= summary_data_sum
  network_summary_dict["Data_P2P_Baseline"] = summary_data

  # calculate communication rates for different parts of the network
  network_rates = [[] for x in xrange(len(network_labels))]
  bytes_data = data_dict["Data_P2P_Bytes"]
  time_data = data_dict["Data_P2P_SendTime"]
  for i in range(nranks):
    for j in range(nranks):
      if bytes_data[i,j] > 0:
        network_rates[int(network_data[i,j])-1].append(bytes_data[i,j]/time_data[i,j])

  for i, label in enumerate(network_labels):
    np_network_rates = np.array(network_rates[i])
    rate_mean = np.mean(np_network_rates)
    rate_std = np.std(np_network_rates)
    print label, " comms rate is ", round(rate_mean,3), " +/- ", round(rate_std,3), " GB/s"


    
def point_size_update(new_point_size):
  global data_all
  global pc_dict

  for key in pc_dict:
    pc_dict[key].set_sizes(np.full(data_all['cnt'], new_point_size))


def totuple(lst):
  tp = ()
  for item in lst:
    tp += (item,)
  return tp
    



if not have_package_pylab():
  print "Unable to display graph as matplotlib could not be imported."
  sys.exit()
  
script_title = "app2mosaic_full_display"
script_version = "v1.0.0"
data_label_prefix = "Data_P2P_"

ifn = ""
idata = 0
node2node = False
internodeonly = False
nodesummary = False
display_output = ""
user_pt_size = 0.0
parse_args()

nranks = 0
nnodes = 0
data_dict = {}
pe_node_xyz = {}
pe_node_id = {}
node_ids = []
parse_pat_data()
calculate_network_distances()
node_data_dict = {}

network_labels = ["on-node", "backplane", "group", "optical"]
network_colors = ["#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"]

network_summary_dict = {}
collate_node_data()
 

print "Parsed %s..." %(ifn)
key_list = []
for key in data_dict:
  key_list.append(get_simple_data_label(key))
key_list.sort()
data_key = ""
for i, key in enumerate(key_list):
  if i == idata:
    data_key = key
  print "  %d: %s" %(i, data_label_prefix+key)
if 0 == len(data_key):
  print "Invalid data index, must be value between 0-", str(len(key_list)-1), "."
else:
  print "Loading ", data_key, " data..."

prepare_plot_data(data_key)
if nodesummary:
  prepare_summary_data()


fig = plt.figure()
fig.canvas.set_window_title(script_title + ' ' + script_version)


if nodesummary:
  ind = np.arange(len(network_summary_dict))
  cum_values = np.zeros(len(network_summary_dict))
  width = 0.55
  halfwidth = width / 2

  bars = ()
  bar_labels = ()
  for i, label in enumerate(network_labels):
    values = []
    for key in sorted(network_summary_dict):
      values.append(network_summary_dict[key][i]*100.0)
    
    bars += (plt.bar(ind, totuple(values), width, bottom=totuple(cum_values), color=network_colors[i]),)
    bar_labels += (label,)

    for j, val in enumerate(values):
      cum_values[j] += val
    
  plt.ylabel('Percentage Breakdown')
  xlabels = []
  for key in sorted(network_summary_dict):
    xlabels.append(get_snappy_data_label(key))
  plt.xticks(ind + halfwidth, totuple(xlabels))

  plt.xlim(ind[0]-halfwidth,ind[-1]+width+halfwidth+1.8)
  plt.ylim(0,100)

  bars = bars[::-1]
  bar_labels = bar_labels[::-1]
  plt.legend(bars, bar_labels)
  plt.tight_layout()

else:
  ax = fig.add_subplot(1,1,1)
  ax.set_aspect('equal', adjustable='box')
  plt.subplots_adjust(left=0.3,bottom=0.3)

  min_pt_size = 10.0 if node2node else 1.0
  max_pt_size = 100.0 if node2node else 20.0
  pt_size = 20.0 if node2node else 1.0
  if user_pt_size > 0.0:
    pt_size = user_pt_size
  pc_dict = {}
  # todo: points that coincide with axis limits should not be truncated
  if "network" in data_key.lower():
    if not node2node:
      pt_size = 2.0
    for i in range(len(network_labels)):
      if internodeonly and network_labels[i] == "on-node":
        continue
      index_sel = np.where(data_all["vals"] == i+1)[0]
      if len(index_sel) > 0:
        pc = plt.scatter(data_all["xpts"][index_sel], data_all["ypts"][index_sel],
               c=network_colors[i], lw=0, s=pt_size, marker='o')
        pc_dict[network_labels[i]] = pc
    recs = []
    for i in range(len(network_colors)):
      recs.append(mpatches.Rectangle((0,0),1,1,fc=network_colors[-1-i]))
    display_network_labels = network_labels[:0:-1] if internodeonly else network_labels[::-1]
    plt.legend(recs, display_network_labels, bbox_to_anchor=(-0.1, 1.03))
  else:
    cm = plt.cm.get_cmap('plasma_r')
    pc = plt.scatter(data_all["xpts"], data_all["ypts"], c=data_all["vals"],
           vmin=data_min, vmax=data_max, lw=0, s=pt_size, marker='o', cmap=cm)
    pc_dict['all'] = pc
    plt.colorbar(pc)


  # setup axes limits, ticks and labels
  nitems = nnodes if node2node else nranks
  xlabel = "target "
  ylabel = "origin "
  xlabel = "node" if node2node else "core"
  ylabel = "node" if node2node else "core"
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position("top")
  ax.set_xlim(0,nitems)
  ax.set_ylim(nitems,0)
  ax.xaxis.set_label_text(xlabel)
  ax.yaxis.set_label_text(ylabel)

  # setup data label
  # todo: centre data labels
  xoff = 0.1 if node2node else 0.1
  yoff = 0.08 if node2node else 0.08
  data_label_xpos = nitems/2 - nitems*xoff
  data_label_ypos = nitems + nitems*yoff
  snappy_label = get_snappy_data_label(data_label)
  if len(data_unit_label) > 0:
    snappy_label += " [" + data_unit_label + "]"
  plt.text(data_label_xpos, data_label_ypos, snappy_label)

  if user_pt_size == 0.0:
    # add slider to control scatter point size
    # todo: improve alignment of point size slider
    pointSizeAxes = plt.axes([0.375, 0.12, 0.45, 0.03])
    pointSizeSlider = Slider(pointSizeAxes, "point size", min_pt_size, max_pt_size, valinit=pt_size, color="black")
    pointSizeSlider.on_changed(point_size_update)

  
plt.show()
if len(display_output) > 0:
  fig.savefig(display_output + ".png", bbox_inches='tight')


