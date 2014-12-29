#!/usr/bin/ruby
$VERBOSE=nil

$timerun='timerun'
$feature='featuresSAT12'
$solver_p1=Array.new()
$solver_p1[0]='' # start with index 1 as in matlab
$solver_p1[1]='ebglucose.sh $CNF'
$solver_p1[2]='ebminisat_static $CNF -verb=0'
$solver_p1[3]='glucose.sh $CNF'
$solver_p1[4]='glueminisat_static -compe $CNF'
$solver_p1[5]='lingeling-static $CNF'
$solver_p1[6]='lr_gl_shr_static $CNF -verb=0 -model'
$solver_p1[7]='minisatpsm.sh $CNF'
$solver_p1[8]='MPhaseSAT64-static  $CNF'
$solver_p1[9]='precosat $CNF'
$solver_p1[10]='qutersat_release $CNF'
$solver_p1[11]='rcl_static.sh $CNF'
$solver_p1[12]='RestartSAT -rfirst=1 -var-decay=0.95 $CNF'
$solver_p1[13]='cryptominisat-snt-st $CNF'
$solver_p1[14]='Spear-32_121 --nosplash --model-stdout --spset-sw-verif --dimacs $CNF --tmout 7200 --seed $SEED'
$solver_p1[15]='Spear-32_121 --nosplash --model-stdout --spset-hw-bmc --dimacs $CNF --tmout 7200 --seed $SEED'
$solver_p1[16]='EagleUP-static $CNF $SEED'
$solver_p1[17]='sparrow2011-static $CNF $SEED'
$solver_p1[18]='march_rw-static $CNF'
$solver_p1[19]='MPhaseSAT_M-static $CNF'
$solver_p1[20]='sattime2011-static $CNF $SEED'
$solver_p1[21]='TNM $CNF $SEED'
$solver_p1[22]='mxc-sat09 -i $CNF'
$solver_p1[23]='gnovelty+2 $CNF $SEED'
$solver_p1[24]='sattime-static $CNF $SEED'
$solver_p1[25]='sattime+-static $CNF $SEED'
$solver_p1[26]='clasp2 --sat-p=20,25,120,-1,0 --berk-max=512 --berk-once --otfs=1 --recursive-str --dinit=800,10000 --dsched=20000,1.1 --reverse-arcs=2  $CNF'
$solver_p1[27]='clasp1 --dimacs --number=1 --sat-p=20,25,150 --hParam=0,512 --file=$CNF'
$solver_p1[28]='picosat-static $CNF'
$solver_p1[29]='MPhaseSAT-static $CNF'
$solver_p1[30]='SApperloT2010-static $CNF'
$solver_p1[31]='sol.sh $CNF'


cnfname=ARGV[0]
seed=ARGV[1]

for i in 1...$solver_p1.length()
   mycmdfoo1 = $solver_p1[i].sub(/\$CNF/, "#{cnfname}") 
    mycmdfoo = mycmdfoo1.sub(/\$SEED/, "#{seed}")
    
   mycmd = "timerun 5 #{mycmdfoo}"
   system mycmd
end

