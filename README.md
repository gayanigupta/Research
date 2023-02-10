# Research

In order to run the k-core code. 
1. Clone the repo. 
2. Use the follwoing commands from the terminal. 
(i) Move to the folder in where your Make File is located. 
/Users/<username>/Documents/GitHub/Research/graph_kcores/graph_kcores/Example-v2.0
(ii) Ensure your Make File has a correct entry. 
kcore: 	test_kcores.cpp      
	g++ -g -Wall $(BASIC_IO_INC) $(BASIC_SETOP_INC) $(BASIC_TRAVERSAL_INC) $(BASIC_CHANGE_INC) $(BASIC_ANALYSIS_INC) test_kcores.cpp -o core
(iii) Execure the command 
make kcore
(iv)next type ./core

You must be able to see the outoput as : 


--Network with only K-cores--
Core 0's vertexs: 
Core 1's vertexs: 
Core 2's vertexs: 3, 4, 6
Core 3's vertexs: 4, 6, 7
Core 4's vertexs: 6, 7
Core 5's vertexs: 
Core 6's vertexs: 7


