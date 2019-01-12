#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;



if ($opt_j == 1) {goto JUMPMARK; }

my $chdPidNb=0;
## == number of proc ==
my $num_proc = 0;
## == number of collected ==
my $num_collect = 0;
my $collect;
## == get the child signal ==    
$SIG{CHLD} = \&proNBminus;
sub proNBminus{ 
	$num_proc--; 
	#print "\naaaaa -- \$num_proc=$num_proc\n"; warn "\naaaaa -- \$num_proc=$num_proc\n";
} 

my $runningPidHash;

 
my $blastOutIdx=1; my $BlastOutFileHash; 
MARKONE: foreach my $algae (    sort { $a cmp $b } (   keys  (  %{ $AllAgaeSeqPathHash }  )   )    ){
	if ( -d ("$outDir/$algae") ){}else{ system ("mkdir -p $outDir/$algae"); }
	foreach my $GenOrEst (    sort { $a cmp $b } (   keys  (  %{ $AllAgaeSeqPathHash->{$algae} }  )   )    ){
		
		foreach my $PROorRNA (    sort { $a cmp $b } (   keys  (  %{ $CisElementsPathHash }  )   )    ){
			
			foreach my $CisElemnt (    sort { $a cmp $b } (   keys  (  %{ $CisElementsPathHash->{$PROorRNA} }  )   )    ){
				if ( -d ("$outDir/$algae/$CisElemnt/$GenOrEst") ){}else{ system ("mkdir -p $outDir/$algae/$CisElemnt/$GenOrEst"); }
				
				
				my $algaeCisDIr="$outDir/$algae/$CisElemnt/$GenOrEst";
				my $algaeCisDIr_link="$outDirBaseNM/$algae/$CisElemnt/$GenOrEst";
				
				my $blastProg='tblastn'; if ($PROorRNA eq 'RNA'){ $blastProg='blastn'; }
				
				my $blastOutPath="$outDir/$algae/$CisElemnt/$GenOrEst/blo.out";  my $baseBlastOutPath="$outDirBaseNM/$algae/$CisElemnt/$GenOrEst/blo.out";  $blastOutIdx++; 
				#if ($blastOutIdx>$testNb) {last MARKONE;} ############################  ²âÊÔÓÃ  ############################ 
				#pipe(PRT_RDR_1, CHD_WTR_1);
				
				my $bestBlastFile="$outDir/$algae/$CisElemnt/$GenOrEst/best.blo.html"; 
				my $BlastOutHashFile="$outDir/$algae/$CisElemnt/$GenOrEst/blo.hash";  #blo.out   blo.hash  
				my $Not_reducedBlOutHashFile="$outDir/$algae/$CisElemnt/$GenOrEst/longBlo.hash";
				
				#$RealhyLinkDir="${opt_o}.HpLkDIR"; my $hyLinkDir
				my $RealCpOutHpLinkBestOutFile="$RealhyLinkDir/${algae}__${CisElemnt}__${GenOrEst}.best.blo.html";
				my $CpOutHpLinkBestOutFile="$hyLinkDir/${algae}__${CisElemnt}__${GenOrEst}.best.blo.html";
				
				$chdPidNb++;
				my $pid=fork();    
				if (!defined($pid)) {
          print "\naaaaa Error in fork: $!";
        exit 1;
        }
        if ($pid == 0) {
          ## == child proc ==
          #print "\naaaaa Child $chdPidNb : My pid = $$\n";  
          warn "\nChild $chdPidNb : My pid = $$\n";
          
          if ($opt_b == 1){
					  warn "blastall -i $CisElementsPathHash->{$PROorRNA}->{$CisElemnt} -F F -a 5 -d $AllAgaeSeqPathHash->{$algae}->{$GenOrEst} -p $blastProg -e $eYuzhi -o $blastOutPath\n";
				    system ("blastall -i $CisElementsPathHash->{$PROorRNA}->{$CisElemnt} -F F -a 5  -d $AllAgaeSeqPathHash->{$algae}->{$GenOrEst} -p $blastProg -e $eYuzhi -o $blastOutPath");
			    }
			    
			    my ( $Not_reducedBlOutHash, $blastOutHash)=@{ &BioPerlBlastPraser_bl2seq_20170330( $AllAgaeSeqPathHash, $CisElementsPathHash, $GenOrEst, $algae, $PROorRNA, $CisElemnt, $algaeCisDIr, $algaeCisDIr_link, $blastOutPath, $bestBlastFile) };  
			    $blastOutHash->{'realPathFile'}=$blastOutPath;   
			    $blastOutHash->{'LinkPathFile'}=$baseBlastOutPath; 
			    $blastOutHash->{'PraseHashFile'}=$BlastOutHashFile;
			    $blastOutHash->{'BestEvResult'}=$bestBlastFile;
			    $blastOutHash->{'BestBlaRstFi'}=$RealCpOutHpLinkBestOutFile;
			    $blastOutHash->{'HpLkBestBlaRstFi'}=$CpOutHpLinkBestOutFile;
			    
			    if ($blastOutHash->{'totalHit'}>0){
			      warn "cp -rf $bestBlastFile $RealCpOutHpLinkBestOutFile\n";
			      system ("cp -rf $bestBlastFile $RealCpOutHpLinkBestOutFile");
			    }
				  my $blastDtBs=$blastOutHash->{'_ResultArray'}->[0]->{'_database_name'}; $blastDtBs=~s/^\s*//g; $blastDtBs=~s/\s*$//g;
				  my $blastInDB=$AllAgaeSeqPathHash->{$algae}->{$GenOrEst};               $blastInDB=~s/^\s*//g; $blastInDB=~s/\s*$//g;
				  my $fistQuryHead=&GetHead                   ( $blastOutHash->{'_ResultArray'}->[0]->{'_query_name'} );  
				  my $inputWord   =&GetFirstLine              ( $CisElementsPathHash->{$PROorRNA}->{$CisElemnt}       );

          if ( (-s ($CisElementsPathHash->{$PROorRNA}->{$CisElemnt}) )>0 ){
					  if ( $blastDtBs eq $blastInDB ){ #if ($blastOutHash->{'_ResultArray'}->[0]->{'_database_name'} eq $AllAgaeSeqPathHash->{$algae}->{$GenOrEst}){
					  	if ($fistQuryHead eq $inputWord){
					  		
					  	  &PrintDumper ( $BlastOutHashFile,        $blastOutHash );
					  	  &PrintDumper ( $Not_reducedBlOutHashFile,$Not_reducedBlOutHash );
				        
				        
			        }
			        else{ print  "\n\n\$fistQuryHead=$fistQuryHead !eq! $inputWord=\$inputWord\n\n";  
			        	die "\n\n\$fistQuryHead=$fistQuryHead !eq! $inputWord=\$inputWord\n\n";
			        }
				      
					  }
					  else {print 	 "$blastDtBs=$blastDtBs eq $blastInDB=\$blastInDB \n\n\$blastOutHash->{'_ResultArray'}->[0]->{'_database_name'}=$blastOutHash->{'_ResultArray'}->[0]->{'_database_name'}\n\$AllAgaeSeqPathHash->{$algae}->{$GenOrEst}=$AllAgaeSeqPathHash->{$algae}->{$GenOrEst}\n";
						  die "\$blastDtBs=$blastDtBs eq $blastInDB=\$blastInDB \n\n\$blastOutHash->{'_ResultArray'}->[0]->{'_database_name'}=$blastOutHash->{'_ResultArray'}->[0]->{'_database_name'}\n\$AllAgaeSeqPathHash->{$algae}->{$GenOrEst}=$AllAgaeSeqPathHash->{$algae}->{$GenOrEst}\n";
					  }
				  }
				  else {
				  	my $tpSize=(  -s ( $CisElementsPathHash->{$PROorRNA}->{$CisElemnt} )  ); 
				  	warn "\n\n\n(-s (\$CisElementsPathHash->{$PROorRNA}->{$CisElemnt}=$CisElementsPathHash->{$PROorRNA}->{$CisElemnt})=$tpSize\n\n\n\n";
				  }

			    
          exit 0;
        }
				 
				warn "\naaaaaa \$pid=$pid\n";   #print "\naaaaaa \$pid=$pid\n"; 
				
				$runningPidHash->{$pid}=1;
				Marktwo: while(1) {
				  my $runningPidNB=0; 
				  my $updatePidHash;
				  foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
				  	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); }#  print "\naaaaaa\$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\n";}
				  	else {$updatePidHash->{$pidKey}=1; $runningPidNB++;                         }# print "\naaaaaa\$runningPidNB child pid is runing!!\t$pidKey=$pidKey is still running!!! \n\n";}
				  }
				  warn "\naaaaaa \$runningPidNB=$runningPidNB\n";   #print "\naaaaaa Totally \$runningPidNB=$runningPidNB is running!\n"; 
				  $runningPidHash=$updatePidHash;
				  if ($runningPidNB<$coreNb2Use){warn"aaaaa \$pid=$pid\t\$runningPidNB=$runningPidNB <  $coreNb2Use !! Countinue to next fork!\n\n";last Marktwo;    }
				  else                          {warn"aaaaa \$pid=$pid\t\$runningPidNB=$runningPidNB >= $coreNb2Use !! Sleep 1 seconds!\n\n"; sleep(1);   }
				}
				


			}
		}
	}
}
my $time_1=`date`;
warn "\n\nStart to wait !!\n\n", $time_1, "\n\n" ;

Markthree: while(1) {
  my $runningPidNB=0; 
  my $updatePidHash;
  foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
  	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); my $ps=`ps $pidKey`; print "\naaaaaa\$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\t\$ps=$ps\n\n"; warn "\naaaaaa\$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\$ps=$ps\n\n";}
  	else {
  		my $ps=`ps $pidKey`; 
  		if ($ps=~m/$pidKey/){
  		  $updatePidHash->{$pidKey}=1; $runningPidNB++;  
  		  print "\naaaaaa $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!! \n\$ps=$ps\n\n"; 
  		  warn "\naaaaaa $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!! \t\$ps=$ps\n\n";
  		}else {warn "\naaaaaa Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps\n\n"; print "\naaaaaa Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps\n\n";}
  	}
  }
  
  warn "\naaaaaa Totally  \$runningPidNB=$runningPidNB\n";   print "\naaaaaa Totally \$runningPidNB=$runningPidNB is running!\n"; 
  $runningPidHash=$updatePidHash;
  if ($runningPidNB>0){ print"aaaaa \$runningPidNB=$runningPidNB > 0 !! Sleep 1 seconds!\tCountinue to next Wait checking run!\n\n"; sleep(1);   }   
  else                 {print"aaaaa \$runningPidNB=$runningPidNB <= 0 !! all child pid dead\n\n";  sleep(1);last Markthree;   }
}

my $waitPid=wait(); my $time_2=`date`; warn "Start Wait AT:\n$time_1\n\nEnd of Wait, \$waitPid=$waitPid!!\n\n", $time_2, "\n\n" ;



#do {            sleep(2); warn "\n\nSleep 2 seconds! \$num_proc=$num_proc\n\n\n";  print  "\n\n\naaaaa Sleep 2 seconds! \$num_proc=$num_proc\n\n\n";        } until ($num_proc < 1);
JUMPMARK: my $allCisEmBlsWiseHash;
my $BlastExcelOutHash;
my @allAgae=@{ &getDirArray($outDir) };
foreach my $eachAg (@allAgae){  my $cisDir="$outDir/$eachAg"; my @allCisElms=@{ &getDirArray($cisDir) };
	foreach my $eachCisElm (@allCisElms){ my $GeEsDir="$outDir/$eachAg/$eachCisElm"; my @allGeEs=@{ &getDirArray($GeEsDir) };
		foreach my $eachGeEs (@allGeEs){
			my $bltFile="$outDir/$eachAg/$eachCisElm/$eachGeEs/blo.out";    if (-e($bltFile)){}  else {die "cannot open $bltFile : $!\n\n";}
			my $hashFile="$outDir/$eachAg/$eachCisElm/$eachGeEs/blo.hash";  if (-e($hashFile)){} else {die "cannot open $hashFile : $!\n\n";}
			my $hashHere;			$hashHere=&retrieve($hashFile);
			
			#my $rdsHashFile="$outDir/$eachAg/$eachCisElm/$eachGeEs/rdsbl.hash";  if (-e($rdsHashFile)){} else {die "cannot open $rdsHashFile : $!\n\n";}
			#my $rdsHashHere;			$rdsHashHere=&retrieve($rdsHashFile);         
			
			#my $rdsHashFile="$outDir/$eachAg/$eachCisElm/$eachGeEs/rdsbl.hash";  if (-e($rdsHashFile)){} else {die "cannot open $rdsHashFile : $!\n\n";}
			#my $hashHere=&retrieve($rdsHashFile);   
			  
			$allCisEmBlsWiseHash->{$eachCisElm}->{$eachAg}->{$eachGeEs}=$hashHere;
			#$allCisEmBlsWiseHash->{$eachCisElm}->{$eachAg}->{$eachGeEs}=$rdsHashHere;
			
			my $bestHit='NotFound';
		  if ($hashHere->{'totalHit'}>0){ #print "\cl\cl\n\n\&print_all_sub_array(\$hashHere)\n\n"; &print_all_sub_array($hashHere); print "\cl\cl\n\n\n\n"; 
		    my $bestQuery=$hashHere->{'_ResultArray'}
		                               ->[ $hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'QurIdx'} ]
		                               ->{'_query_name'}; 
		    my $ShortBestQuery=$bestQuery;   $ShortBestQuery=~s/^(\S+)\s*.*$/$1/;                         
		    my $bestSbj=$hashHere->{'_ResultArray'}
		                             ->[ $hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'QurIdx'} ]				                           
		                             ->{'_hitArray'}
		                             ->[ $hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'HitIdx'} ]
		                             ->{'_accessionNB'};  
		    my $ShortBestSbj=$bestSbj;   $ShortBestSbj=~s/^(\S+)\s*.*$/$1/;                          
		    my $bestEv=$hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'SotVal'};
		    my $bestCoserv=$hashHere->{'_ResultArray'}
		                                ->[ $hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'QurIdx'} ]				                           
		                                ->{'_hitArray'}
		                                ->[ $hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'HitIdx'} ]
		                                ->{'_hitTotalConserved'};  $bestCoserv=$bestCoserv*100; $bestCoserv=sprintf "%.2f",$bestCoserv; $bestCoserv="$bestCoserv%";
		    my $bestIdenty=$hashHere->{'_ResultArray'}
		                                ->[ $hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'QurIdx'} ]				                           
		                                ->{'_hitArray'}
		                                ->[ $hashHere->{'SortHash'}->{'EvalueSort'}->[0]->{'HitIdx'} ]
		                                ->{'_hitTotalIdentical'};  $bestIdenty=$bestIdenty*100; $bestIdenty=sprintf "%.2f",$bestIdenty; $bestIdenty="$bestIdenty%";                          
		    $bestHit="$bestEv($bestCoserv|$bestIdenty)$ShortBestQuery|$bestSbj"; 
		    
		    
		  
		                            
		  }
		  #print "\cl\cl\$BlastOutFileHash->{$algae}->{$GenOrEst}->{$PROorRNA}->{$CisElemnt}->{'RealPath'}=\$blastOutPath=$BlastOutFileHash->{$algae}->{$GenOrEst}->{$PROorRNA}->{$CisElemnt}->{'RealPath'}\n\n\&print_all_sub_array(\$hashHere)\n\n"; &print_all_sub_array($hashHere); print "\cl\cl\n\n\n\n"; 
      if ($bestHit eq 'NotFound'){
      	$BlastExcelOutHash->{$eachAg}->{$eachCisElm}->{$eachGeEs}=$bestHit;
      }
      else{
		    $BlastExcelOutHash->{$eachAg}->{$eachCisElm}->{$eachGeEs}=&makeHyperLk( $hashHere->{'HpLkBestBlaRstFi'}, $bestHit ) ;  
			}
			
		}
	}
}
#$outDir/$algae/$CisElemnt/$GenOrEst
