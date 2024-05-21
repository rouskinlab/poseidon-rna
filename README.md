Functions

	poseidonrna.graph(filename, num_strucs = 1, prefix = '#', color = None, color_type = None, bases = False, node_size = 10)

		filename: str with .ct or .db filename
  		num_strucs: int number of structures to pull from file or string 'all' to graph all
    		prefix: str key used in ct file that is common in the name of the structure to parse rows with new structure title. Not necessary for db file
      		color: str varna coloring filename.
		color_type: str if 'varna' will color using varnafile above
  		bases: boolean if True labels nucleotides. 
    		node_size: how big the nodes are drawn
    		

		output: plot of 2D RNA structure organized by Kamada-Kawai Algorithm
