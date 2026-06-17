cd dynamics
swig -python dynamicsUtils.i
python setup.py build_ext --inplace

cd ..
cd kinematics
swig -python kinematicsUtils.i
python setup.py build_ext --inplace
pause