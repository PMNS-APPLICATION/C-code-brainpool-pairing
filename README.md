# C-code-brainpool-pairing
This repository contains the codes of the PMNS generated for the moduli of the Brainpool and Pairing curves.
<br>
<br>
Each directory contains two main files:
    <br>
  - **main.c**: it allows to check that operations are performed correctly.
    <br>
    Compilation command:
    > gcc -Wall -O3 main.c -o main -lgmp
    
    Execution command:
    > ./main 
    
  - **main__nb_cycles.c**: for timing comparisons (in number of clock cycles).
    <br>
    Compilation command:
    > gcc -Wall -O3 main__nb_cycles.c -o cmain -lgmp -lcrypto
    
    Execution command:
    > ./cmain 
    >

    Note that you need to disable the turbo boost to get accurate timing results.
