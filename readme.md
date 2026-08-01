Project with scripts written to develop a multiscale arterial wall model.

# to do a simulation with abaqus.

1. create the abaqus input file:
  >$ abaqus cae noGUI=test_2d_cae.py

2. run abaqus simulation:
  >$ python3 run_abq.py

3. process abaqus outputadatabase (results) to a text file:
  >$ abaqus python test_odb.py
  
4. create graphics or plots from the test file resulted from the odb reading:
  >$ python3 test_plot.py
  
## this is the testing case, test the UMAT before running complex simulations

### the test is performed with one element simulation to understand possible mistakes

#### ABAQUS user-routines user guide

http://abaqusdocs.eait.uq.edu.au/v6.11/books/sub/default.htm