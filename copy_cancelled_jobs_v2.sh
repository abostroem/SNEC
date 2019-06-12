#USEAGE: bash copy_cancelled_jobs.sh 226 225 224
# awk '{print $10,$1}' cancelled_jobs.txt | sort | cut -d '_' -f 2
# awk '{print "$10,$1"}' cancelled_jobs.txt | sort > copy_cancelled_job_list.txt

export SNECPATH='/home/bostroem/2018an/SNEC-master'
for line in $(awk '{printf("%s+%s\n",$10,$1)}' cancelled_jobs.txt)
do
    job=$(echo "$line" | cut -d '_' -f 2)
    cpu=$(echo "$line" | cut -d '+' -f 1)
    #export SEED=$(awk 'NR==$job {print $0}' input_dir_list.txt)
    export SEED=$(cat $SNECPATH/input_dir_list.txt | head -n $job | tail -n 1)
    echo "copying $job on $cpu to $SEED"
	#echo $SEED
    #Check that this does the right thing
    scp -r $cpu:/scratch/bostroem/job_$job/* $SEED
done
