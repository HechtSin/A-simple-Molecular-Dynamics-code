### A simple molecular dynamics code ##
### Xing He    ###
### PHY 566    ###
### 04/17/2018 ###


######## Running the scripts ################
1. To run NVT by leapfrog

python3 leapfrog_andersen.py

This will generate energy.txt and temerature.txt
To plot them, run:

python3 plot_temperature.py
python3 plot_energy.py

2. To run NVE by leapfrog

comment out Andersen part in leapfrog.py
Then run:

python3 leapfrog_andersen.py

3. To run verlet NVE

python3 verlet.py


####### File descriptions #######

1. initialize.py
To initial atom velocities

2. LJ.py
To define the Lennard-Johns potential
return both potential energy and force

3. leapfrog.py
To run leafrog MD

4. verlet.py
To run verlet MD

5. energy.txt & temperature.txt
Output from MD

6. plot_temperature.py & plot_energy.py
To plot output of MD
