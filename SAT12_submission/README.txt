How to use it
================================
ruby SATzilla12.rb <type> <CNF>
     type: RAND (random SAT category); HAND (hard Combinatorial SAT+UNSAT category); INDU (application  SAT+UNSAT category); ALL (mixture of random SAT+UNSAT, crafted, and application).
     CNF: SAT instance name in CNF formula


Notes:
================================
Please check candidate solvers before you run SATzilla. 
1.	"cd bin"
2.	Run "testsolver.rb satexample.cnf 1234 |grep SATISFIABLE |wc -l", you should get 28
3.	Run "testsolver.rb unsatexample.cnf 1234 |grep UNSATISFIABLE |wc -l", you should get 28
4.	Run "testsolver.rb easy.cnf 1234 |grep SATISFIABLE |wc -l", you should get 31


If you have any question, please contact Lin Xu at xulin730@cs.ubc.ca


Best

Lin

