#!/usr/bin/perl -w
use strict;
use Cwd;
use LWP::Simple;
use FindBin qw($Bin);

my $currdir = getcwd;

my $VERSION_NUM = "2.2.7MR";

require("$Bin\/_getmarker_prodigal.pl");

my $MAX_MARKER = 0;
my $LOGOUT;

my $SETTING_FILE = "setting";
my $HMMSEARCH = "hmmsearch";

checkProgram();

my $RSCRIPT = "Rscript";
my $HEATMAP_R = "$Bin\/heatmap.r";
my $MARKERHMM = "$Bin/marker.hmm";
my $MARKERNUM = 107;
my $MAXBIN = "$Bin\/src\/MaxBin";
my $ABUND_OUTPUT = "";
my $THREADNUM = 1;
my $REASSEMBLY = 0;
my $RPLOT = 0;
my $KMERLEN = 55;
my $RECURSION_MAX = 5;
my $MIN_SEQ_LENGTH = 1000;
my $MIN_BIN_SIZE = 100000;
my $FINISH = "FINISH";

my $USAGE = qq(MaxBin - a metagenomics binning software [Edited by: Malte Ruehlemann (m.ruehlemann\@ikmb.uni-kiel.de)].
Usage:
  run_MaxBin.pl
    -contig (contig file)
    -prodigal (prodigal amino-acid output file)
    -out (output file)

   (Input reads and abundance information)
    [-abund (abundance file) -abund2 (abundfile) -abund3 (abundfile) -abund4 ... ]

   (You can also input lists consisting of reads and abundance files)
    [-abund_list (list of abundance files)]

   (Other parameters)
    [-min_contig_length (minimum contig length. Default 1000)]
    [-max_iteration (maximum Expectation-Maximization algorithm iteration number. Default 50)]
    [-thread (thread num; default 1)]
    [-prob_threshold (probability threshold for EM final classification. Default 0.9)]
    [-plotmarker]
    [-markerset (marker gene sets, 107 (default), 40 (Bacteria & Archaea),
      or 120 (recommended; GTDB Bacteria Marker Genes), and 122 (GTDB Archaea Marker Genes)  See README for more information.)]

  (for debug purpose)
    [-version] [-v] (print version number)
    [-verbose]
    [-preserve_intermediate]

  Please specify -abund information.
  You can input multiple abundance files at the same time.
  Please read README file for more details.\n);


main();

sub main
{
	my $contig_f = "";
	my $contig_name = "";
  my $prodigal_f = "";
  my $prodigal_name = "";
	my @abund_f;
	my $abundcount = 0;
	my $abund_list = "";
	my $out_f = "";
	my $outname = "";
	my $outdir = "";
	my $verbose = 0;
	my $preserve = 0;
	my $maxem = -1;
	my $remodel_contig = 0;
	my $prob_threshold = -1;
  my $old_contig_name;
  my $old_prodigal_name;
	my $new_contig;
  my $new_prodigal;
	my $isfastq;
	my @finisharr;

	my $starttime = time();
	my $endtime;

	my $ARGC = @ARGV;
	my $i;
	my $j;
	my $k;
	my @arr;
	my @tmparr;
	my $cmd;
	my $line;
	my $tmp;
	my $param;
	my $param2;

	# Create a temporary log file
	$k = 0;
	while ($k == 0)
	{
		$k = int(rand(100000000));
		if (-e "$k.log")
		{
			$k = 0;
		}
	}
	open(TMPLOG, ">$k.log");


	# Check if MaxBin program exist. If not then build everything, including 3rd-party software
	if (!(-e $MAXBIN))
	{
		print "Program not built. Please run \"make\" under src directory to build MaxBin program.\n";
		close(TMPLOG);
		unlink("$k.log");
		exit(-1);
	}

	print "MaxBin $VERSION_NUM\n";
	print TMPLOG "MaxBin $VERSION_NUM\n";
	for ($i = 0; $i < $ARGC; $i++)
	{
		if ($ARGV[$i] eq "-contig")
		{
			$i++;
			$contig_f = $ARGV[$i];
			$contig_name = $ARGV[$i];
			print "Input contig: $contig_f\n";
			print TMPLOG "Input contig: $contig_f\n";
		}
    elsif ($ARGV[$i] eq "-prodigal")
    {
      $i++;
      $prodigal_f = $ARGV[$i];
      $prodigal_name = $ARGV[$i];
      print "Input prodigal: $prodigal_f\n";
      print TMPLOG "Input prodigal: $prodigal_f\n";
    }
		elsif ($ARGV[$i] eq "-abund_list")
		{
			$i++;
			$abund_list = $ARGV[$i];
		}
		elsif ($ARGV[$i] =~ /^\-abund/)
		{
			$i++;
			if (-e $ARGV[$i])
			{
				$abund_f[$abundcount] = $ARGV[$i];
				$abundcount++;
				print "Located abundance file [$ARGV[$i]]\n";
				print TMPLOG "Located abundance file [$ARGV[$i]]\n";
			}
			else
			{
				print "Cannot find abundance file [$ARGV[$i]]\n";
				close(TMPLOG);
				unlink("$k.log");
				exit;
			}
		}
		elsif ($ARGV[$i] eq "-out")
		{
			$i++;
			$out_f = $ARGV[$i];
			$j = rindex($out_f, "/");
			if ($j == -1)
			{
				$outname = $out_f;
			}
			else
			{
				$outname = substr($out_f, $j + 1);
				$outdir = substr($out_f, 0, $j);
			}
			print "out header: $ARGV[$i]\n";
			print TMPLOG "out header: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-thread")
		{
			$i++;
			$THREADNUM = $ARGV[$i];
			print "Thread: $ARGV[$i]\n";
			print TMPLOG "Thread: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-max_iteration")
		{
			$i++;
			$maxem = $ARGV[$i];
			print "Max iteration: $ARGV[$i]\n";
			print TMPLOG "Max iteration: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-prob_threshold")
		{
			$i++;
			if ($ARGV[$i] >= 0)
			{
				$prob_threshold = $ARGV[$i];
				print "Probability threshold: $ARGV[$i]\n";
				print TMPLOG "Probability threshold: $ARGV[$i]\n";
			}
			else
			{
				print "Invalid value $ARGV[$i]. Probability threshold not changed.\n";
			}
		}
		elsif ($ARGV[$i] eq "-min_contig_length")
		{
			$i++;
			$MIN_SEQ_LENGTH = $ARGV[$i];
			print "Min contig length: $ARGV[$i]\n";
			print TMPLOG "Min contig length: $ARGV[$i]\n";
		}
		# elsif ($ARGV[$i] eq "-reassembly") # Reassembly currently still does not support FASTQ reads
		# {
		# 	$REASSEMBLY = 1;
		# 	print "Reassembly: 1\n";
		# 	print TMPLOG "Reassembly: 1\n";
		# }
		elsif ($ARGV[$i] eq "-markerset")
		{
			$i++;
			if ($ARGV[$i] == 40)
			{
				print "Switch to 40 marker genes universal for bacteria and archaea.\n";
				print TMPLOG "Switch to 40 marker genes universal for bacteria and archaea.\n";
				$MARKERHMM = "$Bin/bacar_marker.hmm";
				$MARKERNUM = 40;
			}
      ### Adding GTDB bac120 and ar122 sets
      if ($ARGV[$i] == 120)
			{
				print "Switch to 120 marker genes universal for bacteria from the GTDB.\n";
				print TMPLOG "Switch to 120 marker genes universal for bacteria from the GTDB.\n";
				$MARKERHMM = "$Bin/gtdb_bac120.hmm";
				$MARKERNUM = 120;
			}
      if ($ARGV[$i] == 122)
      {
        print "Switch to 122 marker genes universal for archaea from the GTDB.\n";
        print TMPLOG "Switch to 122 marker genes universal for archaea from the GTDB.\n";
        $MARKERHMM = "$Bin/gtdb_ar122.hmm";
        $MARKERNUM = 122;
      }

		}
		elsif ($ARGV[$i] eq "-verbose")
		{
			$verbose = 1;
		}
		# elsif ($ARGV[$i] eq "-plotmarker")
		# {
		# 	$RPLOT = 1;
		# }
		elsif ($ARGV[$i] eq "-preserve_intermediate")
		{
			$preserve = 1;
		}
		elsif ($ARGV[$i] eq "-version" || $ARGV[$i] eq "-v")
		{
			print "MaxBin $VERSION_NUM\n";
			exit;
		}
		else
		{
			print "Unrecognized token \[$ARGV[$i]\]\n";
			print $USAGE;
			close(TMPLOG);
			unlink("$k.log");
			exit;
		}
	}

	if ($contig_f eq "")
	{
		print "No Contig file. Please specify contig file by -contig\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}
	if ($out_f eq "")
	{
		print "Please specify output file by -out.\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}
	if (-e "$out_f.log")
	{
		unlink "$out_f.log";
	}

	if ($abund_list ne "")
	{
		open(FILE, "<$abund_list") || die "Cannot open abundance list file $abund_list\n";
		while(defined($line = <FILE>))
		{
			chomp($line);
			if (-e $line)
			{
				print "Located abundance file [$line]\n";
				print TMPLOG "Located abundance file [$line]\n";
				$abund_f[$abundcount] = $line;
				$abundcount++;
			}
			elsif ($line ne "")
			{
				print "Cannot find abundance file [$line]. Stop.\n";
				exit;
			}
		}
		close(FILE);
	}

	if ($abundcount == 0)
	{
		print "Please input at least one abundance file. You may also specify an abundance file list.\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}

	$new_contig = checkContig($contig_f, $MIN_SEQ_LENGTH, "$out_f.contig.tmp", "$out_f.tooshort");
	if ($new_contig ne "")
	{
		$old_contig_name = $contig_f;
		$contig_f = $new_contig;
		$remodel_contig = 1;
	}

  $new_prodigal = checkProdigalFile($prodigal_f, "$out_f.prodigal.tmp");
  if ($new_prodigal ne "")
  {
    $old_prodigal_name = $prodigal_f;
    $prodigal_f = $new_prodigal;
  }

	$param = "";
	$j = 1;
	for ($i = 0; $i < $abundcount; $i++)
	{
		$line = "$contig_f.abund" . $j;
		checkAbundFile($abund_f[$i], $line);
		if ($j == 1)
		{
			$param = $param . " -abund " . $line;
		}
		else
		{
			$param = $param . " -abund$j " . $line;
		}
		$j++;
	}

	openLOG("$out_f.log");
	close(TMPLOG);
	open(TMPLOG, "<$k.log");
	while(defined($line = <TMPLOG>))
	{
		print $LOGOUT $line;
	}
	close(TMPLOG);
	unlink("$k.log");


	my @binarr = ();
	my $currbin;
	my $currout;
	my $currnum;
	my $maxnum;
	my @summarr = ();
	my @allsumarr = ();
	my @reclassifyarr = ();
	my @noclassarr = ();
	push(@binarr, $contig_f);
	# Push some dummy into result array for the original dataset
	push(@summarr, "dummy");
	push(@allsumarr, "dummy");
	push(@reclassifyarr, 1);
	$currnum = 0;
	$maxnum = 1;

	# Run HMMER3 to identify seed contigs
	writeLOG("Searching against $MARKERNUM marker genes to find starting seed contigs for [$prodigal_f]...\n");
	if (!(-e "$prodigal_f.hmmout.$FINISH"))
	{
		#if (-e "$prodigal_f")
		#{
		#	unlink("$prodigal_f");
		#}
		getHMM($prodigal_f, "$prodigal_f.hmmout");
		touch("$prodigal_f.hmmout.$FINISH");
	}

	while ($currnum < $maxnum)
	{
		$currbin = $binarr[$currnum];

		if ($currbin eq $contig_f)
		{
			$currout = $out_f;
		}
		else
		{
			$currout = substr($currbin, 0, length($currbin) - 6) . ".out";
		}
		if ($MAX_MARKER == 1)
		{
			$i = gethmmmarker("$prodigal_f.hmmout", $currbin, $MIN_SEQ_LENGTH, "$currout.seed", $MARKERNUM, 1);
		}
		else
		{
			$i = gethmmmarker("$prodigal_f.hmmout", $currbin, $MIN_SEQ_LENGTH, "$currout.seed", $MARKERNUM);
		}
		# Check if current file exceeds the current file processing recursion limitation
		# Check for how many "####.out", in which # represent digits
		@tmparr = $currbin  =~ /[0-9]{4}.out/g;
		$j = scalar @tmparr;

		if ($currbin eq $contig_f && $i == -1)
		{
			writeLOG("Try harder to dig out marker genes from contigs.\n");
			$i = gethmmmarker("$prodigal_f.hmmout", $currbin, $MIN_SEQ_LENGTH, "$currout.seed", $MARKERNUM, 1);
			if ($i == -1)
			{
				writeLOG("Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop.\n");
				closeLOG("$out_f.log");
				exit(-1);
			}
		}
		elsif ($currbin ne $contig_f && $j >= $RECURSION_MAX)
		{
			$currnum++;
			next;
		}
		elsif ($currbin ne $contig_f && $i == -1)
		{
			$currnum++;
			next;
		}
		#elsif ($j > $RECURSION_MAX)
		#{
		#	$currnum++;
		#	next;
		#}
		elsif ($i != -1)
		{
			$reclassifyarr[$currnum] = 1;
		}

		# Running MaxBin
		writeLOG("Done data collection. Running MaxBin...\n");
		$param2 = "";
		if ($maxem != -1)
		{
			$param2 = "-max_run $maxem ";
		}
		if ($prob_threshold != -1)
		{
			$param2 = $param2 . "-prob_threshold $prob_threshold ";
		}
		if ($verbose == 1)
		{
			$param2 = $param2 . "-verbose ";
		}
		if ($THREADNUM > 1)
		{
			$param2 = $param2 . "-thread $THREADNUM";
		}
		$cmd = "$MAXBIN -fasta $currbin $param -seed $currout.seed -out $currout -min_contig_length $MIN_SEQ_LENGTH $param2";
		writeLOG("Command: $cmd\n");
		if (!(-e "$currout.$FINISH"))
		{
			system($cmd);
		}

		if (checkResult("$currout.summary") == -1)
		{
			if ($currbin eq $contig_f)
			{
				writeLOG("Error encountered while running core MaxBin program. Error recorded in $currout.log.\nProgram Stop.\n");
				exit(-1);
			}
		}
		push(@noclassarr, $currout);

		# Read in summary file and putting bins into stack
		open(FILE, "<$currout.summary");
		while(defined($line = <FILE>))
		{
			chomp($line);
			if ($line =~ /^Bin \[([A-Za-z0-9._\\\/\@\!\|\#\$\%\^\?\<\>\[\]\{\}\(\)\+\-]+)\]\t([0-9.\t]+)/)
			{
				push(@binarr, $1);
				@arr = split(/\t/, $2);
				push(@summarr, $arr[0]);
				push(@allsumarr, $2);
				push(@reclassifyarr, 0);
				$maxnum++;
			}
		}
		close(FILE);

		$currnum++;
		touch("$currout.$FINISH");
		push(@finisharr, "$currout.$FINISH");
	}

	# Remove progress recording files that ends with $FINISH

	unlink("$prodigal_f.hmmout.$FINISH");
	foreach $i (@finisharr)
	{
		unlink($i);
	}

	# Put everything in order

	# Remove re-classified bins
	my @tmpgenomesize = ();
	my @tmpgc = ();
	my @genomesize = ();
	my @gc = ();
	my @bin_to_noclass = ();
	for ($k = 0; $k < $maxnum; $k++)
	{
		($tmpgenomesize[$k], $tmpgc[$k]) = getBinInfo($binarr[$k]);
	}
	for ($k = 0; $k < $maxnum; $k++)
	{
		if ($reclassifyarr[$k] == 1 || $tmpgenomesize[$k] < $MIN_BIN_SIZE)
		{
			if ($reclassifyarr[$k] == 1 && $binarr[$k] ne $contig_f)
			{
				unlink($binarr[$k]);
			}
			elsif ($reclassifyarr[$k] == 0 && $tmpgenomesize[$k] < $MIN_BIN_SIZE)
			{
				push(@bin_to_noclass, $binarr[$k]);
			}
			splice(@reclassifyarr, $k, 1);
			splice(@binarr, $k, 1);
			splice(@allsumarr, $k, 1);
			splice(@summarr, $k, 1);
			splice(@tmpgenomesize, $k, 1);
			splice(@tmpgc, $k, 1);
			$k--;
			$maxnum--;
		}
	}

	# First sort @summarr and store the index in another array
	my @sortarr = sort {$b <=> $a} @summarr;

	#my %indexhash;
	#@indexhash{@sortarr} = (0..$#summarr);

	my @sorted;
	my @indexarr;

	for ($i = 0; $i < $maxnum; $i++)
	{
		for ($j = 0; $j < $maxnum; $j++)
		{
			if (!(defined($sorted[$j])) && $summarr[$i] == $sortarr[$j])
			{
				$indexarr[$i] = $j;
				$sorted[$j] = 1;
				last;
			}
		}
	}
	@sortarr = ();

	my @filearr;
	for ($k = 0; $k < $maxnum; $k++)
	{
		#$i = $indexhash{$summarr[$k]};
		$i = $indexarr[$k];
		$j = getBinNum($i + 1, $maxnum);
		# rename original fasta file
		rename ($binarr[$k], "$out_f.$j.fasta.tmp");
		$filearr[$i] = "$out_f.$j.fasta.tmp";
		$genomesize[$i] = $tmpgenomesize[$k];
		$gc[$i] = $tmpgc[$k];
		$sortarr[$i] = $allsumarr[$k];
	}
	# Write summary
	open(TMPSUM, ">$out_f.tmp.summary");
	for ($i = 0; $i < $maxnum; $i++)
	{
		$j = getBinNum($i + 1, $maxnum);
		print TMPSUM "Bin \[$out_f.$j.fasta\]\t$sortarr[$i]\n";
	}
	close(TMPSUM);

	# Collect all no-classes and log file
	open(OUT, ">$out_f.tmp.noclass");
	$j = @noclassarr;
	for ($i = 0; $i < $j; $i++)
	{
		open(FILE, "<$noclassarr[$i].noclass");
		while(defined($line = <FILE>))
		{
			print OUT $line;
		}
		close(FILE);

		open(FILE, "<$noclassarr[$i].log");
		while(defined($line = <FILE>))
		{
			writeLOG($line, 0);
		}
		writeLOG("\n");
		close(FILE);
	}
	foreach $currbin (@bin_to_noclass)
	{
		if (-e $currbin)
		{
			open(FILE, "<$currbin");
			while(defined($line = <FILE>))
			{
				print OUT $line;
			}
			close(FILE);
			unlink($currbin);
		}
		else
		{
			writeLOG("File $currbin not found.\n");
		}
	}
	close(OUT);

	# delete all normal files
	#for ($i = 0; $i < $maxnum; $i++)
	#{
	#	if ($binarr[$i] ne $contig_f && (-e $binarr[$i]))
	#	{
	#		unlink($binarr[$i]);
	#	}
	#}
	$j = @noclassarr;
	for ($i = 0; $i < $j; $i++)
	{
		unlink("$noclassarr[$i].log");
		unlink("$noclassarr[$i].summary");
		unlink("$noclassarr[$i].noclass");
		#unlink("$noclassarr[$i].prob_dist");
		#unlink("$noclassarr[$i].dist");
		unlink("$noclassarr[$i].seed");
	}

	# rename all tmp files to normal files
	rename("$out_f.tmp.summary", "$out_f.summary");
	rename("$out_f.tmp.noclass", "$out_f.noclass");
	for ($i = 0; $i < $maxnum; $i++)
	{
		$currbin = substr($filearr[$i], 0, length($filearr[$i]) - 4);
		rename($filearr[$i], $currbin);
		$filearr[$i] = $currbin;
	}


	# Count marker genes for bins
	my $completearr;

	$completearr = countmarker("$prodigal_f.hmmout", $out_f, $outname, $MARKERHMM, "$out_f.marker", $MARKERNUM);
	$i = 0;

	# Re-read summary file and write other info: completeness, Genome size, GC
	open(FILE, "<$out_f.summary");
	open(SUMOUT, ">$out_f.summary.tmp");
	$j = $abundcount;
	if ($j == 1)
	{
		print SUMOUT "Bin name\tAbundance\tCompleteness\tGenome size\tGC content\n";
	}
	else
	{
		print SUMOUT "Bin name\tCompleteness\tGenome size\tGC content\n";
		open(SUMABUND, ">$out_f.abundance");
		print SUMABUND "Bin name";
		for ($i = 0; $i < $abundcount; $i++)
		{
			print SUMABUND "\t$abund_f[$i]";
		}
		print SUMABUND "\n";
	}
	$i = 0;
	while(defined($line = <FILE>))
	{
		if ($line =~ /^Bin \[$out_f.([0-9]+).fasta\]\t([0-9.\t]+)/)
		{
			if ($j == 1)
			{
				printf SUMOUT "$outname.$1.fasta\t%0.2f\t%0.1f%%\t$genomesize[$i]\t%0.1f\n", $2, $$completearr[$i] * 100, $gc[$i];
			}
			else
			{
				print SUMOUT "$outname.$1.fasta";
				printf SUMOUT "\t%0.1f%%\t$genomesize[$i]\t%0.1f\n", $$completearr[$i] * 100, $gc[$i];
				print SUMABUND "$outname.$1.fasta";
				@arr = split(/\t/, $2);
				foreach $tmp (@arr)
				{
					printf SUMABUND "\t%0.2f", $tmp;
				}
				print SUMABUND "\n";
			}
			$i++;
		}
	}
	close(FILE);
	close(SUMOUT);
	if ($j > 1)
	{
		close(SUMABUND);
	}
	rename("$out_f.summary.tmp", "$out_f.summary");

	if ($RPLOT == 1)
	{
		$cmd = "$RSCRIPT $HEATMAP_R $out_f.marker $out_f.marker.pdf";
		print "$cmd\n";
		system($cmd);
	}

	# Get markers for each bin
	listmarker_bybin("$prodigal_f.hmmout", $out_f, "$prodigal_f", $MARKERHMM, $out_f, $MARKERNUM);
	if ($outdir ne "")
	{
		chdir($outdir);
	}
	$cmd = "tar zcvf $outname.marker_of_each_bin.tar.gz $outname.*.marker.fasta";
	system($cmd);
	$cmd = "rm $outname.*.marker.fasta";
	system($cmd);
	if ($outdir ne "")
	{
		chdir($currdir);
	}

	if ($preserve == 0)
	{
		writeLOG("Deleting intermediate files.\n");

		$cmd = "rm -rf $out_f.tmp";
		system($cmd);
		#unlink("$contig_f.prob_dist");
		#unlink("$contig_f.seed.m8");
		unlink("$contig_f.seed");
		for ($i = 0; $i < $abundcount; $i++)
		{
			$j = $i + 1;
			$line = "$contig_f.abund" . $j;
			unlink($line);
		}
		unlink("$prodigal_f.hmmout");
	}

	if ($remodel_contig == 1)
	{
		unlink($contig_f);
	}

	# Write post instruction to users
	writeLOG("\n\n========== Job finished ==========\nYielded $maxnum bins for contig (scaffold) file $contig_name\n\n");
	writeLOG("Here are the output files for this run.\nPlease refer to the README file for further details.\n\n");
	$i = $maxnum;
	$j = getBinNum($i, $maxnum);
	$tmp = getBinNum(1, $maxnum);
	writeLOG("Summary file: $out_f.summary\n");
	$k =  $abundcount;
	if ($k > 1)
	{
		writeLOG("Genome abundance info file: $out_f.abundance\n");
	}
	writeLOG("Marker counts: $out_f.marker\nMarker genes for each bin: $out_f.marker_of_each_gene.tar.gz\nBin files: $out_f.$tmp.fasta - $out_f.$j.fasta\nUnbinned sequences: $out_f.noclass\n");

	if ($RPLOT == 1 && -e "$out_f.marker.pdf")
	{
		writeLOG("Marker plot: $out_f.marker.pdf\n");
	}
	writeLOG("\n\n========== Elapsed Time ==========\n");
	$endtime = time();
	$line = getElapsedTime($endtime - $starttime);
	writeLOG("$line\n");

	closeLOG("$out_f.log");
}


sub getHMM
{
	my $prodigal_f = $_[0];
	my $out_f = $_[1];
	my $cutmethod = $_[2];
	my $cmd;

	print "Running HMMER hmmsearch....\n";
	if ($MARKERNUM == 40)
	{
    $cmd = "$HMMSEARCH --domtblout $out_f -E 1e-3 --cpu $THREADNUM $MARKERHMM $prodigal_f 1>$out_f.out 2>$out_f.err";
	}
	else
	{
    $cmd = "$HMMSEARCH --domtblout $out_f --cut_tc --cpu $THREADNUM $MARKERHMM $prodigal_f 1>$out_f.out 2>$out_f.err";
	}
	system($cmd);
	if (checkResult("$out_f") == -1)
	{
		print "Error running Hmmer3. Output recorded in $out_f.out and $out_f.err.\nProgram Stop.\n";
		exit(-1);
	}
	else
	{
		unlink("$out_f.out");
		unlink("$out_f.err");
	}

}

sub getBinNum
{
	my $snum = $_[0];
	my $smax = $_[1];
	my $tmp = "";
	my $ret = "";
	if ($smax > 1000)
	{
		if ($snum < 10)
		{
			$tmp = "000" . $snum;
		}
		elsif ($snum < 100)
		{
			$tmp = "00" . $snum;
		}
		elsif ($snum < 1000)
		{
			$tmp = "0" . $snum;
		}
		else
		{
			$tmp = $snum;
		}
	}
	else
	{
		if ($snum < 10)
		{
			$tmp = "00" . $snum;
		}
		elsif ($snum < 100)
		{
			$tmp = "0" . $snum;
		}
		else
		{
			$tmp = $snum;
		}
	}

	#if ($smax > 1000 && $snum < 1000)
	#{
	#	$ret = "0" . $tmp;
	#}
	#else
	#{
		$ret = $tmp;
	#}
	return $ret;
}

sub checkResult
{
	my $s = -s $_[0];
	if (-e $_[0] && $s > 0)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

sub checkContig
{
	my $contig_f = $_[0];
	my $min_length = $_[1];
	my $faout = $_[2];
	my $below_f = $_[3];
	my $i;
	my $len;
	my $faline;
	my $header;
	my $seq;
	my $gz = "";
	my $FA_LIMIT = 128;
	my $FASTA_LEN = 70;

	if ($contig_f =~ /.gz$/)
	{
		$gz = $contig_f . ".fa";
		$i = "gunzip -c $contig_f > $gz";
		system($i);
		$contig_f = $gz;
	}

	open(FILE, "<$contig_f") || die "Cannot open file $contig_f\n";
	$i = 0;

	# Remodel input contig file
	seek(FILE, 0, 0);
	$header = "";
	#$faout = $contig_f . ".tmp";
	open(FAOUT, ">$faout") || die "Cannot write into specified output file/directory. Please check your settings or disk space.\n";
	open(FABELOW, ">$below_f");
	while(defined($faline = <FILE>))
	{
		chomp($faline);
		if ($faline =~ /^>/)
		{
			if ($header ne "" && ($len = length($seq)) >= $min_length)
			{
				print FAOUT "$header\n";
				$i = 0;
				while ($i < $len)
				{
					if ($i + $FASTA_LEN > $len)
					{
						print FAOUT substr($seq, $i);
						print FAOUT "\n";
						$i = $len;
					}
					else
					{
						print FAOUT substr($seq, $i, $FASTA_LEN);
						print FAOUT "\n";
						$i += $FASTA_LEN;
					}
				}
			}
			elsif ($header ne "" && ($len = length($seq)) < $min_length)
			{
				print FABELOW "$header\n";
				$i = 0;
				while ($i < $len)
				{
					if ($i + $FASTA_LEN > $len)
					{
						print FABELOW substr($seq, $i);
						print FABELOW "\n";
						$i = $len;
					}
					else
					{
						print FABELOW substr($seq, $i, $FASTA_LEN);
						print FABELOW "\n";
						$i += $FASTA_LEN;
					}
				}
			}
			$header = $faline;
			$seq = "";
		}
		else
		{
			$seq = $seq . $faline;
		}
	}
	$len = length($seq);
	if ($len >= $min_length)
	{
		print FAOUT "$header\n";
		$i = 0;
		while ($i < $len)
		{
			if ($i + $FASTA_LEN > $len)
			{
				print FAOUT substr($seq, $i);
				print FAOUT "\n";
				$i = $len;
			}
			else
			{
				print FAOUT substr($seq, $i, $FASTA_LEN);
				print FAOUT "\n";
				$i += $FASTA_LEN;
			}
		}
	}
	else
	{
		print FABELOW "$header\n";
		$i = 0;
		while ($i < $len)
		{
			if ($i + $FASTA_LEN > $len)
			{
				print FABELOW substr($seq, $i);
				print FABELOW "\n";
				$i = $len;
			}
			else
			{
				print FABELOW substr($seq, $i, $FASTA_LEN);
				print FABELOW "\n";
				$i += $FASTA_LEN;
			}
		}
	}
	close(FILE);
	close(FAOUT);
	close(FABELOW);

	if ($gz ne "")
	{
		unlink($gz);
	}

	return $faout;
}

# (Genome size, GC content) = getBinInfo(fasta file)
sub getBinInfo
{
	my $all = 0;
	my $gc = 0;
	my $i;
	my $line;
	open(FAFILE, "<$_[0]");
	while(defined($line = <FAFILE>))
	{
		if ($line !~ /^>/)
		{
			chomp($line);
			$i = ($line =~ tr/[cgCG]//);
			$gc += $i;
			$i = ($line =~ tr/[atcgATCG]//);
			$all += $i;
		}
	}
	close(FAFILE);

	if ($all == 0)
	{
		$gc = 0;
	}
	else
	{
		$gc = $gc / $all * 100;
	}
	return ($all, $gc);
}

# (Read and Check program directory)
sub checkProgram
{
	my $line;
	my $tmpstr;
	my $tmpname =  "tmp_" . time();

	open(FILE, "<$Bin\/$SETTING_FILE");
	while(defined($line = <FILE>))
	{
		chomp($line);
		if ($line =~ /\[([A-Za-z0-9_]+)\] ([A-Za-z0-9._\(\)\[\]\{\}\|\$\!\=\-\+\\\/]+)/)
		{
			if ($1 eq "HMMER3")
			{
				if (-d $2 && -e "$2\/$HMMSEARCH")
				{
					$HMMSEARCH = $2 . "\/" . $HMMSEARCH;
				}
			}

		}
	}
	close(FILE);

	# HMMER3
	$line = "$HMMSEARCH 1>$tmpname 2>/dev/null";
	system($line);
	$tmpstr = "";
	open(FILE, "<$tmpname");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /hmmsearch/)
	{
		print "Cannot run HMMER3. Please indicate the file directory in \'setting\' file.\n";
		exit;
	}

	unlink("$tmpname");
}

sub checkProdigalFile
{
	my $prodigal_f = $_[0];
	my $prodigal_out_f = $_[1];
	my $line;
  my @splitline;

	if ($prodigal_f eq $prodigal_out_f)
	{
		return;
	}

	#Check header
	open(PRODIN, "<$prodigal_f") || die "Cannot open abund file $prodigal_f\n";
	open(PRODOUT, ">$prodigal_out_f") || die "Cannot open tmp abund file $prodigal_out_f\n";
	while(defined($line = <PRODIN>))
  {
		if ($line =~ /^>/)
		{
      @splitline = split(' ', $line);
			print PRODOUT "$splitline[0]\n";
		}
		else
		{
			print PRODOUT $line;
		}
	}
	close(ABUNDOUT);
	close(PRODIN);
  return($prodigal_out_f)
}

sub checkAbundFile
{
	my $input_f = $_[0];
	my $out_f = $_[1];
	my $line;

	if ($input_f eq $out_f)
	{
		return;
	}

	#Check header
	open(ABUNDIN, "<$input_f") || die "Cannot open abund file $input_f\n";
	open(ABUNDOUT, ">$out_f") || die "Cannot open tmp abund file $out_f\n";
	while(defined($line = <ABUNDIN>))
	{
		if ($line =~ /^>/)
		{
			print ABUNDOUT substr($line, 1);
		}
		else
		{
			print ABUNDOUT $line;
		}
	}
	close(ABUNDOUT);
	close(ABUNDIN);
}



sub getElapsedTime
{
	my $input = $_[0];
	my $hour;
	my $min;
	my $sec;
	my $str;

	$sec = $input % 60;
	$input = int($input / 60);
	$min = $input % 60;
	$input = int($input / 60);
	$hour = $input;

	$str = $hour . " hours " . $min . " minutes and " . $sec . " seconds.\n";
	return $str;
}

sub openLOG
{
	my $log_f = $_[0];
	open($LOGOUT, ">$log_f.tmp");
}

# Currently the writeLOG procedure will NOT write to STDOUT only if there is a second parameter, say 1. Will be revised later.
sub writeLOG
{
	my $msg = $_[0];
	my $is_print = $_[1];
	if (!(defined($is_print)))
	{
		print $msg;
	}
	print $LOGOUT $msg;
}

sub closeLOG
{
	my $log_f = $_[0];
	close($LOGOUT);
	rename("$log_f.tmp", $log_f);
}

sub touch
{
	my $cmd = "touch $_[0]";
	system($cmd);
}
