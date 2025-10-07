These is a SLIM modified version of the MWIS code 

https://www.stha.de/shares/mwis2025/mwis2025_code.tar.zst,

its version used for METAMIS compariosn,  from 
 
 mwis2024_code/integer/variant4-dualgap-heuristicspeedup-mindual=0.1-integer/ 

 folder. Here we

 - modified  include/mpopt/mwis/solver_cont_temp.hpp -> solver_cont_temp.hpp
 - added main.cpp to read json files
 - simple CMAkeLists.txt for compilation

1. Compile the MWIS solver: 

  1.1. Download json library from https://github.com/nlohmann/json/releases and put it into /external. 
       Check the respective path to <nlohmann/json.hpp> in CMakeLIsts.txt

  1.2. Initialize cmake: 

		  mkdir build && cd build
		  cmake ..
		  cmake -DCMAKE_BUILD_TYPE=Release ..
		  cmake --build .
   

2. Run the solver by calling:
   
   ./main <input-file.json> <number-of-greedy-soluitons>  > <output-file-name>

   2.1 To unpack json.xz files use

      xz -d input-file.json.xz

   2.2  Examples:  

   ./main test3.json 16 > output3.txt
   ./main AM_030.json 10 > output4.txt
  

	2.3 Note: 

	1) The second parameter defines how many greedy solutions will be generated. The program stops after producing these greedy solutions.
	2) The precision of solving the relaxed problem is defined in main.cpp in the line with the comment "relative duality gap" of

	    solver.run(  50,          // batch size
                1000000,      // max number of batches
                0.1           // relative duality gap
    );


	Change this value (0.1) for a trade-off between time and solutions quality.

	3) The code for generating and printing the greedy solutions is now in main.cpp. So, now you do not need to print it, just generate and copy into a given container with mwis::solver.generate_greedy_solution(begin,end). 
	Note that the generated solution is represented in a form of a binary vector.

3. Run the extract_solutions.py script to extract a json-file from the output:

	python3 extract_solutions.py <output-file-name> <json-file-name>

	Examples: 

	python3 extract_solutions.py output3.txt result3.json
	python3 extract_solutions.py output4.txt result4.json