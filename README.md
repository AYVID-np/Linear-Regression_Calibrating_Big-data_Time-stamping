# Linear-Regression_Calibrating_Big-data_Time-stamping
Calibrating the experimentally derived "big data"  to generate a "time-stamped" dataset employing "linear Regression Fit" on 792 degrees of freedom.

-> The following script is written in C using CERN ROOT libraries.
-> This project involved calibrating the big data (~ TB) derived out of an experiment into a "time-stamped" data. The time-stamping of the data is in nanoseconds which is the nuclear reaction timescale.
-> The peaks in a 12-Bit Time-to-Digital converter are identified with a resolution of 0.2%. Further, a Linear Regression model is used to fit 792 experimental degrees of freedom.
-> The extracted fitted parameters are employed in calibrating the measurements and assigning a time-stamping to the measured nuclear reaction products
for around e^12 events.
