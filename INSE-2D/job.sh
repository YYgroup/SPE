echo "#!/bin/sh">submit.sh
echo "yhrun -N 1 -n 16 -p work ./DNS >out 2>&1">>submit.sh
chmod 700 submit.sh 
yhbatch -N 1 -n 16 -p work -J SPE-NS_test.mzy  ./submit.sh
rm submit.sh