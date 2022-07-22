#!/usr/bin/perl

#SortData.pl

#This script takes the output from FairShapley.pl and transforms it into a tab separated table. 

#Modules needed for running the script
use strict;
use warnings;


#For each data set, the following variables need to be adjusted:
#	my $results = Combined list of FP indices for gene and species trees obtained from running FairShapley.pl in the following order: 
#	FP indices on gene trees, average FP index across gene trees (possibly manually corrected, see FairShapley.pl), FP index on species tree 
#	my $output = name of the desired output file
#	my @species_names = list of all species names occuring in the data set
#	my $num_trees = number of gene trees used when running FairShapley.pl
#


#Dolphin data set:
#my $results = 'dolphin_fp_combined.txt';
#my $output =  'dolphin_fp_combined_table.txt';
#my @species_names = ("St_bre", "Gr_gri", "Gl_mac", "Pe_ele", "Fe_att", "Ps_cra", "De_cap", "De_del", "St_lon", "Tu_adu", "Tu_tru", "St_coe", "So_chi", "La_hos", "St_att", "So_flu", "Ce_com", "Li_bor", "La_obl", "La_obs", "Or_orc", "La_alb", "La_acu", "Ne_pho", "Ph_pho", "De_leu", "Mo_mon", "Ph_mac");
#my $num_trees = 22;


#Fungi data set:
#my $results = 'fungi_fp_combined.txt';
#my $output =  'fungi_fp_combined_table.txt';
#my @species_names = ("cjadii", "wanosu", "kmarsu", "scerre", "hlib91", "hvalsi", "hsin91", "hjak91", "khat91", "hlac70", "hpse91", "hopuea", "hguiii", "hgui70", "hnec31", "htha91", "hmey91", "hcle91", "huva70", "huva08", "huva86", "huva9-", "han_91", "hocc29", "hocc91", "hosm91", "hosmal", "hvinea", "hvin91");
#my $num_trees = 683;


#Mammal data set:
#my $results = 'mammals_fp_combined.txt';
#my $output =  'mammals_fp_combined_table.txt';
#my @species_names = ("Gal", "Mac", "Mon", "New", "Gor", "Hom", "Pan", "Pon", "Cal", "Tar", "Oto", "Mic", "Tup", "Spe", "Dip", "Mus", "Rat", "Cav", "Ory", "Och", "Sor", "Eri", "Tur", "Bos", "Sus", "Vic", "Fel", "Can", "Equ", "Myo", "Pte", "Lox", "Pro", "Ech", "Das", "Cho", "Orn");
#my $num_trees = 447;


#Plant data set:
#my $results = 'plant_fp_combined.txt';
#my $output = 'plant_fp_combined_table.txt';
#my @species_names = ("aur_pec", "lan_tib", "pet_vol", "pau_tom", "cal_ame", "pro_las", "wes_fru", "col_can", "per_fru", "hyp_sua", "ple_bar", "oci_bas", "lav_ang", "mel_off", "per_atr", "ros_off", "sal_off", "sal_his", "lyc_ame", "men_pip", "men_spi", "mon_did", "ori_maj", "ori_vul", "thy_vul", "pru_vul", "nep_cat", "nep_mus", "gle_hed", "hys_off", "aga_foe", "aju_rep", "teu_can", "cle_bun", "rot_myr", "bal_pse", "mar_vul", "lam_alb", "leo_leo", "leo_car", "phl_fru", "bet_off", "pog_cab", "pet_bam", "hol_san", "scu_bai", "cor_pyr", "con_tom", "vit_agn", "gme_phi", "pre_mic", "tec_gra");
#my $num_trees = 318;


#Primate data set:
#my $results = 'human_chimp_gorilla_fp_combined.txt';
#my $output =  'human_chimp_gorilla_fp_combined_table.txt';
#my @species_names = ("Human", "Chimp", "Gorilla", "Orang");
#my $num_trees = 52;


#Rattlesnake data set:
#my $results = 'rattlesnake_fp_combined.txt';
#my $output =  'rattlesnake_fp_combined_table.txt';
#my @species_names = ("sced127AZ", "sca151WI", "scter115K", "sms1OK", "smc1NC", "smi100FL", "ac1OUTG");
#my $num_trees = 16;


#Rodent data set:
#my $results = 'animal_fp_combined.txt';
#my $output =  'animal_fp_combined_table.txt';
#my @species_names = ("Rat37", "fusci30", "Mus36", "gliro2", "minah5", "major12", "nouhu11", "macro23", "golia7", "ruemm22", "barba33", "lanos14", "roths13", "sp20", "elega10", "myoid34", "eller24", "asper19", "chrys6", "monck4", "argur35", "forre8", "fuscu18", "orali29", "fumeu28", "deser27", "austr26", "fuscu15", "alboc25", "gould17", "condi9", "penic3", "caudi32", "levip21", "saleb31", "rufes16", "imita1");
#my $num_trees = 761;


#Snake data set:
#my $results = 'snake_fp_combined.txt';
#my $output =  'snake_fp_combined_table.txt';
#my @species_names = ("I1_Anolis_carolinensis_XXX", "I33_Python_molurus_XXX", "I7_Acrochordus_granulatus_RAP0485", "I19_Xenodermus_javanicus_FMNH230073", "I12_Pareas_carinatus_CAS247982", "I26_Echis_carinatus_RS142", "I27_Azemiops_feae_KU312227", "I23_Agkistrodon_contortrix_TJL135", "I2_Crotalus_adamanteus_XXX", "I20_Crotalus_horridus_FTB969", "I5_Gerarda_prevostiana_CAS204972", "I3_Lapemis_curtus_CAS204975", "I32_Calliophis_melanurus_RS148", "I4_Lycodryas_inornatus_RAN39035", "I6_Aparallactus_capensis_AMB8180", "I11_Atractaspis_bibroni_AMB8268", "I8_Pseudaspis_cana_AMB6804", "I16_Prosymna_frontalis_MCZ38642", "I15_Psammophis_notostictus_AMB4565", "I24_Hormonotus_modestus_MVZ252503", "I28_Gonionotophis_klingi_MVZ252476", "I10_Sibynophis_bistrigatus_CAS214081", "I25_Grayia_smithii_MCZ182627", "I13_Masticophis_taeniatum_DBS395", "I14_Calamaria_pavimentata_CAS235364", "I17_Pseudoxenodon_macrops_CAS245413", "I9_Diadophis_punctatus_FTB554", "I21_Heterodon_nasicus_MHP8535", "I22_Nerodia_rhombifer_DBS353", "I31_Regina_rigida_FTB594", "I30_Storeria_storerioides_CAS238520", "I18_Storeria_dekayi_FTB2739", "I29_Storeria_occipitomaculata_FTB718");
#my $num_trees = 333;


#Yeast data set:
my $results = 'Yeast_fp_combined.txt';
my $output =  'yeast_fp_combined_table.txt';
my @species_names = ("Sklu","Scer", "Spar","Smik","Skud","Sbay","Scas","Calb");
my $num_trees = 106;



#Open the output file
open(OUTPUT, ">>", $output) or die ("Could not open file $output!\n");

#Print 1st row: Species	Tree1	Tree2	...	Treenumtree	Average	SpeciesTree
print OUTPUT "Species";
for (my $i=1; $i<=$num_trees; $i++){
    print OUTPUT "\tTree".$i;
}
print OUTPUT "\tAverage ";
print OUTPUT "\tSpeciesTree";
print OUTPUT "\n";



#Read in the FP indices for each species and print them to the output file; if a tree T did not contain species x, print NaN for the FP index of x on T
my $tree_counter=0;

foreach my $name (@species_names){
    print OUTPUT $name;

    open(RESULTS, "<", $results) or die ("Could not open file $results!\n");
    while(<RESULTS>){
        my $line = $_;
        chomp $line;
		
		if ($line =~ /tree/){
			$tree_counter++;
		}

        if ($line =~ /$name/ && $tree_counter==1){
            print $line."\n";
            my @tab_split = split(/\t/,$line);
            my $value = $tab_split[1];
            print OUTPUT "\t".$value;
			$tree_counter=0;
        }
		
		elsif ($line =~ /$name/ && $tree_counter>1){
			my @tab_split = split(/\t/,$line);
            my $value = $tab_split[1];
            print $value."\n";
			
			while ($tree_counter > 1){
				print OUTPUT "\t"."NaN";
				$tree_counter--;
			}
			print OUTPUT "\t".$value;
			$tree_counter=0;
		}
        
    }
    
    while ($tree_counter > 0){
        print OUTPUT "\t"."NaN";
        $tree_counter--;
    }
            
    
    print OUTPUT "\n";
    close(RESULTS);
}


close(OUTPUT);
