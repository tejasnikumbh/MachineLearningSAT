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

# ====================================================================
def runSolver(sid, cnfname, cutoff, binpath, inseed)
    seed=rand()
    ssid=sid.to_i
    success=10

    if ssid > $solver_p1.length()-1
       puts "Error, no such solver! #{ssid} > #{$solver_p1.length()}"
       return -1
    end
    tmpfile="/tmp/SATzilla_solver#{seed}"
    mycmdfoo = "#{binpath}/#{$timerun} #{cutoff} #{binpath}/#{$solver_p1[ssid]} > #{tmpfile}"
    mycmdfoo1 = mycmdfoo.sub(/\$CNF/, "#{cnfname}") 
    mycmd = mycmdfoo1.sub(/\$SEED/, "#{inseed}")
    system mycmd
    outflag=0
    if File.exist?(tmpfile)
      checkcmd="grep '^s SATISFIABLE' #{tmpfile} |wc -l"
      File.popen(checkcmd){|lines|
          line=lines.gets
          if line.to_i >0
             outflag=1
          end
       }
      checkcmd="grep '^s UNSATISFIABLE' #{tmpfile} |wc -l"
      File.popen(checkcmd){|lines|
          line=lines.gets
          if line.to_i >0
             puts 's UNSATISFIABLE'
             File.delete(tmpfile)
             return success
          end
       }
      if outflag >0
         File.open(tmpfile){|lines|
           while line=lines.gets
               puts line
           end
         }
         File.delete(tmpfile)
         return success
      end
    end
    return -1
end
# ======================= Main ==========================================
if ARGV.length() <2
   puts "Usage ZILLA.rb <type> <cnf>"
   puts "               -  type: RAND, HAND, INDU, ALL"
   puts "               -  cnf: SAT instance in CNF (must contain p cnf)"
   exit
end
if ARGV.length()==2
 inseed=1234
else
 inseed=ARGV[2].to_i
end

type=ARGV[0]
cnfname=ARGV[1]

execdir=File.dirname(__FILE__)
# the location of bin
binpath="#{execdir}/bin"
modelpath="#{execdir}/models"
featcutoff=90
timefeat=[6, 21, 42, 53, 59, 78, 97, 109, 121, 124];
selected=Array.new()
if type=~/RAND/
 selected=[17, 16, 20, 25,     23, 24] # directly from matlab
 backupm=1    # directly from matlab
 pre1m=1   # directly from matlab
 pre2m=1    # directly from matlab
 pretime1=2.18   # directly from matlab
 pretime2=0    # directly from matlab
end
if type=~/HAND/
 selected=[23, 29, 31, 26, 17, 13, 22, 18, 20, 8, 27, 9]
 backupm=2 # change from 23 to 29 
 pre1m=2
 pre2m=1
 pretime1=4.39
 pretime2=0
end
if type=~/INDU/
 selected=[5, 10, 8, 12, 13, 27, 4, 26, 21, 7, 3]
 backupm=1
 pre1m=7
 pre2m=1
 pretime1=9.86
 pretime2=0
end
if type=~/ALL/
 selected=[5,     17, 26, 31, 13, 18, 20,         3,     8]
 backupm=1
 pre1m=1
 pre2m=1
 pretime1=0
 pretime2=0
end

featmodel="#{type}featureModel.ser"
pwmodel="#{type}pairwiseModel.ser"

ecode=-1
backup=selected[backupm-1] # index start with 0 in ruby
pre1=selected[pre1m-1]
pre2=selected[pre2m-1]
forever=1200 # the longest runtime > 900
# first get the easy feature with 0 cost
mycmd="grep -m 1 '^p cnf' #{cnfname}"
seed=rand()
tmppath='/tmp/'
foo=""
easyfeatname="#{tmppath}easyfeat#{seed}.csv"
featname="#{tmppath}feat#{seed}.csv"
rawfeatname="#{tmppath}rawfeat#{seed}.csv"
# puts easyfeatname
tmpf=File.new(easyfeatname, 'w')
File.popen(mycmd){|lines|
 line=lines.gets
 foo=line.split(" ")
 tmpf.puts "#{foo[2]},#{foo[3]}"
}
tmpf.close

# make prediction if feature is too expensivce
foo=''
mycmd="java -classpath #{modelpath}:#{modelpath}/regression/fastrf SingleClassificationTester #{modelpath}/#{featmodel} #{easyfeatname}"
# puts mycmd
puts "c make prediction of feature cost ..."
featExpen=1
File.popen(mycmd){|lines|
  while line=lines.gets
     foo=line      
#     puts foo
    end
}
featExpen=foo.to_i
File.delete(easyfeatname)

if featExpen==2
  puts "c Feature computation would be too expensive, just run backup solver #{backup}"
  ecode=runSolver(backup, cnfname, forever, binpath, inseed)
     
else
  if pretime1 > 0.1
  puts "c run presolver #{pre1} with cutoff #{pretime1} and check the results"
  ecode=runSolver(pre1, cnfname, pretime1, binpath, inseed)
  if ecode.to_i ==10
     exit
  end
  end
  if pretime2 > 0.1
  puts "c Pre1 fails, run presolver #{pre2} with cutoff #{pretime2} and check the results"
    ecode=runSolver(pre2, cnfname, pretime2, binpath, inseed)
  if ecode.to_i ==10
     exit
  end
  end
  puts "c Presolvers fail, compute features ...."
  
  mycmd="#{binpath}/#{$timerun} #{featcutoff} #{binpath}/#{$feature} -base -sp -dia -cl -ls -lobjois #{cnfname} > #{rawfeatname}"  
  system mycmd
  if File.exist?(rawfeatname)
     mycmd="tail -n 5 #{rawfeatname}"
     # now we need check the features
     foo=[] 
     File.popen(mycmd){|lines|
      while line=lines.gets
          foo=line.split(",")
      end
     }
   #  puts foo.join(",")       
     File.delete(rawfeatname)
     if foo.length()< 125
        puts "c Feature computation did not compute all features, run backup solver #{backup}"
        ecode=runSolver(backup, cnfname, forever, binpath, inseed)
     else
       takefeat=(0..124).to_a - timefeat 
       newfeat=Array.new();
       ind=0;
       for j in 0...takefeat.length()
          foo1= "#{foo[takefeat[j]]}"
       #  no feature should take a negative value; only exception: missing value (-512), treat these as 0  
          newfeat[j]=[foo1.to_f, 0].max
          ind=ind+1
       end    
       puts "c the number of features is #{ind}"
       tmpf=File.new(featname, 'w')
       tmpf.puts newfeat.join(",")
       tmpf.close
       
       puts "c predict which solver should be used ...."
       mycmd="java -classpath #{modelpath}:#{modelpath}/regression/fastrf PairwiseClassificationTester #{modelpath}/#{pwmodel} #{featname}"
 #      puts mycmd
       picksolver=Array.new();
       foo="";
       
       File.popen(mycmd){|lines|
         while line=lines.gets
            if line =~/Rankings\s(.*)/
             foo=$1
             puts "c solver ranking #{foo}"
            end
         end
       }
       if File.exist?(featname)
           File.delete(featname)
       end
       picksolver=foo.split(" ")
         
      # note that selected solver produced by JAVA code start with 0, so we do not need adjust the index
       puts "c run and check best solver #{selected[picksolver[0].to_i]} ..."
     ecode=runSolver(selected[picksolver[0].to_i], cnfname, forever, binpath, inseed)
       if ecode.to_i ==10
        exit
       end
       puts "c run and check second best solver #{selected[picksolver[1].to_i]} ..."
       ecode=runSolver(selected[picksolver[1].to_i], cnfname, forever, binpath, inseed)
       if ecode.to_i ==10
        exit
       end
       
       # in case top 2 solvers fail and we still have time, let's run backup solver
       puts "c top 2 solvers carshed, run backup solver #{backup} ..."
       ecode=runSolver(backup, cnfname, forever, binpath, inseed)
     end
  else
     puts "c featue computation failed, run backup solver #{backup} ..."
     ecode=runSolver(backup, cnfname, forever, binpath, inseed)

  end
end

