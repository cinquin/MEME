//
// simple-shared-doc.js
//
// Function to replace the innerHTML of element "id" with the HTML indicated by "doc_type".
// Easier to read and update than the more flexible approach in shared-doc.js. 
//
function print_doc(id, doc_type) {
  switch (doc_type) { 
    case 'motif-consensus':
      document.getElementById(id).innerHTML = `
	<p> 
	   A consensus sequence is constructed from each column in a
	   motif's frequency matrix using the <b>"50% rule"</b>
	   as follows.
	</p>
	<ol>
	  <li>The letter frequencies in the column are sorted in decreasing order.</li>
	  <li>Letters with frequency less 50% of the maximum are discarded.</li>
	  <li>The letter used in this position in the consensus sequence is determined
	  by the first rule below that applies:</li>
	  <ul>
	    <li>If there is only one letter left, or if the remaining letters exactly match
	    an ambiguous symbol in the alphabet, the <b>letter</b> or <b>ambiguous symbol</b>,
	    respectively, is used.</li>
	    <li>Otherwise, if the remaining set contains at least 50% of the core
	    symbols in the alphabet, the alphabet's <b>wildcard</b>
	    (e.g., "N" for DNA or RNA, and "X" for protein) is used.</li>
	    <li>Otherwise, the letter with the <b>maximum frequency</b> is used.</li>
	  </ul>
	</ol>
      `;
      break;
        
    default:
        document.getElementById(id).innerHTML = "Error--Unrecognized doc_type: " + doc_type;
  }
}
