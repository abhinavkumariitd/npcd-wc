# NPCD-WC
Neighbourhood Proximity based Commmunity Detection using Weighted Centrality (NPCD-WC)
An algorithm designed to detect overlapping community structures in weighted networks.

To be able to use this code, first compile it using the C++ compiler on a Linux/Unix machine as follows:
g++ -std=c++11 npcd_wc.cpp -o npcd_wc

Once the executable file docwan is created, it can be used as follows:

./npcd_wc ./network.txt -w -ov [option] -rh [option]

The file network.txt contains the input network in edge list format. The code accepts three optional flags: -w, -ov, and -rh. The -w flag specifies that the network is weighted; it can be omitted for unweighted networks. The -ov flag sets the overlap threshold (default: 0.4), meaning any two communities with an overlap greater than or equal to this value will be merged. Similarly, the -rh flag sets the neighborhood proximity threshold (default: 0.4). If the user is unsure about the appropriate values for these parameters, they can simply run the code for weighted networks using the default settings.
 If the user is not sure about the values of these parameters, he or she can simply call the code for weighted networks as follows:

./npcd_wc ./network.txt -w            

If the network is unweighted, then the code can be called as follows:

./npcd_wc ./network.txt 

The output is written in the file called npcd_wc-coms.txt, which resides in the same directory where the code lies.


If you use this code for research purpose, please cite the following article:

@article{kumar2024revisiting,
  title={Revisiting neighbourhood proximity based algorithm for overlapping community detection in weighted networks},
  author={Kumar, Abhinav and Kumar, Pawan and Dohare, Ravins},
  journal={Social Network Analysis and Mining},
  volume={14},
  number={1},
  pages={105},
  year={2024},
  publisher={Springer}
}
