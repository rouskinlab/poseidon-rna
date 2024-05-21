import networkx as nx
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
def read_db_structures(filename):
    #only reads first structure
    # Open the .db file
    with open(filename, 'r') as file:
        # Read the lines
        lines = file.readlines()

    names = []
    seqs = []
    pairs = []
    for line in lines:
        if line[0] == '>':
            names.append(line[1:])
        elif line[0] in ['A','C','G','U','a','c','g','u']:
            seq = list(line)[:-1]
        elif line[0] in ['.','(',')']:
            seqs.append(seq)
            pairs.append(list(line)[:-1])

    return names,seqs,pairs

def read_ct_structures(file,prefix):

    # Open the .ctfile
    with open(file, 'r') as file:
        # Read the lines
        lines = file.readlines()

    #Make list of structures with each structure being a list of lines split into the CT columns.
    S= []
    Stemp = ['Null']
    names = []
    for line in lines:

        #separate structures based on some notable "prefix" that appears in each structre title ('#' in Seismic)
        if prefix in line:

            S.append(Stemp)
            Stemp = []
            names.append(line.split()[1])

        else:
            Stemp.append(line.split())

    S.append(Stemp)
    file.close()

    S = S[1:]

    return(names,S)

def create_db_graph(seq,pair,colors):
    #generate empty graph
    G = nx.Graph()

    #create node at infinity
    G.add_node(0,base='N',label='O',color = colors[0])

    #add nodes for each base
    for x in range(len(seq)):

        G.add_node(x+1,base=seq[x],color = colors[x+1])

    #create "pending" list to pair appropriately.
    pending=[]

    for x in range(len(pair)):
        G.add_edge(x,x+1) #from 0 at infinity to 1, then 1 to 2 etc. up to last base

        #Adds a connection point is open brackett
        if pair[x] == '(':
            pending.append(x+1)

        #Connects node to last open brackett from the pending list. 
        elif pair[x] == ')':
            G.add_edge(pending[-1],x+1)
            pending = pending[:-1]

        #ignores if not a brackett i.e. a dot
        else:
            pass

    #G.add_edge(x+1,0)
        
    return G

def conGraphStep(G,line, ROCAUC, color):

    G.add_node(int(line[0]), base = line[1], label = ROCAUC, color = color)

    G.add_edge(int(line[0]),int(line[2]))

    if int(line[0]) < int(line[4]):

        G.add_edge(int(line[0]),int(line[4]))

    return

def create_ct_graph(struc,AUC,colors):
    #generate empty graph
    G = nx.Graph()

    #create node at infinity
    G.add_node(0,base='N',label = 'O',color = colors[0])


    if isinstance(AUC, np.ndarray):# == False:
        for x in range(len(struc)):

            #Adds new node for each base and relevant connected edges
            conGraphStep(G,struc[x],AUC[x],colors[x+1])


    else:
        for x in range(len(struc)):

            #Adds new node for each base and relevant connected edges
            conGraphStep(G,struc[x],'',colors[x+1])

        
    return G

def plotGraph(G,name,color,node_size,bases,view,**kwargs):

    args = dict(kwargs.items())

    plt.rcParams["figure.figsize"] = (15,15)

    try:
        pos = args['pos']
        
    except:
        pos = nx.kamada_kawai_layout(G)

    if color == None:
        
        nx.draw(G,pos,with_labels=False,node_size = node_size)

    else:
        node_colors = [G.nodes[n]['color'] for n in G.nodes]
        nx.draw(G,pos,with_labels=False,node_size = node_size,node_color=node_colors,cmap=sns.color_palette("viridis", as_cmap=True))

    if bases == True:
        #Label nodes with base
        labels = nx.get_node_attributes(G, 'base')
        nx.draw_networkx_labels(G, pos, labels)

    plt.savefig('Structures/'+name+'.svg')
    plt.savefig('Structures/'+name+'.png')

    if view == True:
        plt.show()

    plt.close()

    return(pos)

    
################################################################################

#Create Graph from file
def graph(filename,**kwargs):

    newpath = r'Structures' 
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    args = dict(kwargs.items())

    try:
        num_strucs = args['num_strucs']
        
    except:
        num_strucs = 1

    try:
        prefix = args['prefix']
        
    except:
        prefix = '#'

    #plot color default:none
    try:
        color = args['color']
        color_type = 'varna'
    except:
        color = None
        color_type = None

    #plot bases default:false
    try:
        bases = args['bases']
    except:
        bases = False

    #node size default:10
    try:
        node_size = args['node_size']
    except:
        node_size = 10

    #view graph default:10
    try:
        view = args['view']
    except:
        view = False


    #read in file if dot-bracket
    if filename[-2:]=='db':
    
        names,seqs,pairs = read_db_structures(filename)

    #read in file if connectivity table
    elif filename[-2:] == 'ct':
        
        names, seqs = read_ct_structures(filename,prefix)

    #read in colors from varna file
    if color_type == 'varna':
        with open(color, 'r') as file:
            # Read the lines
            colors = file.readlines()
        newcolors = [0.25]
        for x in range(len(colors)):
            newcolors.append((float(colors[x][:-2])+1)/2)
        file.close()

    else:
        newcolors = np.zeros(len(seqs[0])+1)


    graphs = []
    
    if filename[-2:]=='db':
        
        if num_strucs == 'all':
            num_strucs = len(pairs)
            
        for x in range(num_strucs):
            G = create_db_graph(seqs[x],pairs[x],newcolors)
            graphs.append(G)


    elif filename[-2:] == 'ct':

        if num_strucs == 'all':
            num_strucs = len(seqs)
        
        for x in range(num_strucs):
            G = create_ct_graph(seqs[x],'',newcolors)
            graphs.append(G)

    
    positions = []
    for x in range(num_strucs):

        pos = plotGraph(graphs[x],names[x],color,node_size,bases,view)
        
        positions.append(pos)

    return
