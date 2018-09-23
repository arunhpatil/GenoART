<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>GenoART</title>
<link rel="stylesheet" type="text/css" href="CSS/header2017.css">
<link rel="stylesheet" type="text/css" href="CSS/round_borders.css">
<style type="text/css">
<!--
.style1 {color: #FF0000}
-->
</style>
<style >
#subheader {
    clear: both;
    border-top: 1px dotted #888;
    border-bottom: 1px dotted #888;
    //background: #f8eccf;
    background: #eaeaea;
    //color: #505050;
    padding: 1em;
    margin: 15px 0px 10px 0px;
}
</style>
</head>

<body>
<ul class="headertab">
  <li><a class="active" href="http://GenoART.inhouseprotocols.com/" >Home</a></li>
  <li><a href="http://GenoART.inhouseprotocols.com/faq.php" >FAQ</a></li>
  <li><a href="http://GenoART.inhouseprotocols.com/protocols.php" >Protocols</a></li>
  <li><a href="http://GenoART.inhouseprotocols.com/Contact.php" >Contact</a></li>
  <li><a href="http://GenoART.inhouseprotocols.com/about.php">About</a></li>
  <li style="float:right"><a href="http://ibioinformatics.org/" target="_blank"><em>Institute of Bioinformatics</em></a></li>
</ul>	
<h2 align="center"><em>GenoART: Genome Annotation and Refinement Tool</em></h2>
<table height="69%" width="1298" id="border-radius" class="border-radius-example">

<!--<tr  height="20" valign="top">
   <td height="23" >&nbsp;</td>
   <td colspan="2">&nbsp;</td>
   <td>&nbsp;</td>
 </tr> -->
<tr valign="top">
 
   <td width="3%" height="24" >&nbsp;</td>
   <td colspan="2"> <p><em><strong>Genome Annotation and Refinement Tool</strong></em></p>    <p align="justify" style="font-size: 16px;line-height: 1.30;font-weight: 500;letter-spacing: 0em;font-family: SF Pro Display,SF Pro Icons,Helvetica Neue,Helvetica,Arial,sans-serif;">GenoART is an online tool for detecting coding potential of transcripts complemented by experimentally derived proteomics (MS/MS) data. GenoART database consists of tryptic peptides derived from theoretically translating transcripts into three frames (see Protocols for more). In the current version, a custom database for <em>Homo sapiens</em> is created using transcripts obtained from GENCODE (Release 26, GRCh38.p10) which consists of several coding and non-coding RNA (see table below). This custom database can be used as a reference protein database to search unmatched/unassigned spectra in any search engines such as Mascot, SEQUEST, X! Tandem or Andromeda.  GenoART can be queried to categorize and browse these resultant peptides at the click of a button. The peptides will include the genomic coordinates which allows user to navigate to UCSC genome browser which can facilitate more details about the gene/transcript and can be compaired with predicted genes, Expressed sequence Tag (EST), and perform comparitive genomics for more evidence. A GTF (Gene Transfer Format) will be available to download upon querying, which can later be viewed in any genome visualizer such as UCSC or IGV to compare accross different studies. Download <a href="./Databases/GENCODE.tar"><font color="#0099FF">GENCODE Database</font></a>.</p>
<br />

     <p><em><strong>Preparatory or Preliminary analysis:</strong></em></p> 	
<!--<form action="ihp_proteogenomics_analysis.php"  target="_blank" method="post" enctype="multipart/form-data"> -->
<form action="ihp_proteogenomics_analysis_prelim.php" target="_blank"  method="post" enctype="multipart/form-data"> 
<p align="justify" style="font-size: 16px;line-height: 1.30;font-weight: 500;letter-spacing: 0em;font-family: SF Pro Display,SF Pro Icons,Helvetica Neue,Helvetica,Arial,sans-serif;">In the preparatory analysis, the peptides provided by the user is filtered with the known reference database and remove ambiguous peptides introduced by isobaric ions. A tab delimited text file should be provided as input, where the first column should be peptide sequence. The preparatory analysis yields two files: one mapping to reference proteins and one that doesn't map, and they are named as ProteinCoding and GSSPs respectively. </br></br>NOTE: The file name should not contain special characters or spaces. Please find a sample input file for reference.</p>
          <div id="subheader" align="left">Please select the protein database: 
            <select name="organism">
			<option selected="selected" value="HsRefSeq83.fasta">Human RefSeq 83</option>
		<!--	<option selected="selected" value="./Databases/HsRefSeq83.fasta">Human RefSeq 83</option> -->
		<!--	<option selected="selected" value="HsRefSeq70_cont.fasta">Human RefSeq 70</option> -->
	  <!--              <option value="HsRefSeq65_cont.052614.fasta">Human RefSeq 65</option> -->
            </select>
            <br />
            <br />
         Filename:
            <input type="file" name="file" id="file">
            <input type="submit" name="submit" value="Submit">
            <input type="reset" name = "reset" value="Reset">
            <br />
	</div>
      </form><br />
	      

     <p><em><strong>Analyze genome search-specific peptides:</strong></em></p> 	
<!--<form action="ihp_proteogenomics_analysis.php"  target="_blank" method="post" enctype="multipart/form-data"> -->
<form action="ihp_proteogenomics_analysis.php" target="_blank"  method="post" enctype="multipart/form-data"> 
<p align="justify" style="font-size: 16px;line-height: 1.30;font-weight: 500;letter-spacing: 0em;font-family: SF Pro Display,SF Pro Icons,Helvetica Neue,Helvetica,Arial,sans-serif;">The GSSPs obtained from the preparatory analysis should be provided as input. If you have already performed the preparatory analysis in a different method, you can just load in the list of peptides in this analysis. However, the input should remain as a tab delimited text file, where the first column is peptide sequence. These peptides will be queried against the GenoART database and the results are tabulated in an user friendly browser interface. The results can also be downloaded for further analysis.</br></br>NOTE: The file name should not contain special characters or spaces. Please find a sample input file for reference.</p>
          <div id="subheader" align="left">Upload GSSP file for categorization: 
          <!--  <select name="organism">
			<option selected="selected" value="HsRefSeq83.fasta">Human RefSeq 83</option>
            </select> -->
            <br />
            <br />
         Filename:
            <input type="file" name="file" id="file">
            <input type="submit" name="submit" value="Submit">
            <input type="reset" name = "reset" value="Reset">
            <br />
	</div>
      </form><br />
	 <br />
	 <br />
	      
	</td>
   <td width="4%">&nbsp;</td>

 </tr>
</table>
</body>
</html>
