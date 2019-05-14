#grep ">" | sed -E -e 's/type://g' -e 's/len://g' -e 's/[:()]/ /g' -e 's/-/ /g' | awk '{printf "%s\t%s\tna\t%s,%.2f,%s\t%d \n",$1, $5, $5, sqrt((($8-$7)/($5*3)*100)^2), $4, sqrt(($8-$7)^2)}'


cat cdna.cds | grep ">" | sed -E -e 's/type://g' -e 's/len://g' -e 's/[:()]/ /g' -e 's/-/ /g' | awk '{printf "%s\t%s\tna\t%s,%.2f,%s\t%d \n",$1, $5, $5, sqrt((($8-$7)/($5*3)*100)^2), $4, sqrt(($8-$7)^2)}' > tmp
cat cdna.aa | fasta2oneline.ba | grep -v ">" | sed 's/[^X]//g' | awk '{print length($0)}' > tmp2

paste tmp tmp2 | awk '{print $1, $2, $6, $4, $5}' |column -t 
rm tmp
rm tmp2
