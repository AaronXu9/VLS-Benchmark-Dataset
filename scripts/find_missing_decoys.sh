count=0
missing=0
for dir in /project2/katritch_223/aoxu/PLATE-VS/data/chembl_affinity/uniprot_*; do
    uid=$(basename "$dir" | sed 's/uniprot_//')
    if [ -f "${dir}/${uid}_chembl_activities_filtered.parquet" ]; then
        count=$((count + 1))
        if [ ! -f "${dir}/deepcoy_output/${uid}_generated_decoys.txt" ]; then
            echo "Missing: $uid"
            missing=$((missing + 1))
        fi
    fi
done
echo "Total targets with filtered parquet: $count"
echo "Total missing decoys: $missing"
// Missing: O43598
// Missing: O43613
// Missing: O43614
// Missing: O60259
// Missing: O75417
// Missing: P00492
// Missing: P01024
// Missing: P01106
// Missing: P03951
// Missing: P05067
// Missing: P06858
// Missing: P11362
// Missing: P15309
// Missing: P21802
// Missing: P22455
// Missing: P22607
// Missing: P29597
// Missing: P31639
// Missing: P31751
// Missing: P36776
// Missing: P36897
// Missing: P43490
// Missing: P51449
// Missing: P53779
// Missing: P54764
// Missing: P55072
// Missing: P78540
// Missing: Q01973
// Missing: Q02127
// Missing: Q09472
// Missing: Q13822
// Missing: Q5VWK5
// Missing: Q7Z5W3
// Missing: Q86W56
// Missing: Q92630
// Missing: Q92918
// Missing: Q96GD4
// Missing: Q96P20
// Missing: Q99558
// Missing: Q99683
// Missing: Q9H4A3
// Missing: Q9HD26
// Missing: Q9Y463
// Missing: Q9Y5Z0