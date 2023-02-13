
# BLAST Query Tool

A versatile and user-friendly tool for performing BLAST searches on DNA or protein sequences against multiple databases, including GenBank, RefSeq, and UniProt.

## Features

-   Choose from a variety of BLAST programs (blastn, blastp, blastx, tblastn, or tblastx) for your search.
-   Set parameters such as e-value cutoff, maximum number of hits, and filtering options to customize your search.
-   Get results in a summary report, graphical visualization, and detailed report with information about matches and similarities between sequences.
-   Perform batch BLAST searches on a sequence file in FASTA format.
-   Perform custom BLAST searches with custom parameters.

## Requirements

-   Python 3
-   BioPython

## Usage

The tool provides several functions for performing BLAST searches, generating reports, visualizing results, and exporting the results in different formats.

Here is a brief overview of the functions:

-   `blast_search`: Perform a BLAST search on a single input sequence.
-   `summary_report`: Generate a summary report of the BLAST search results.
-   `detailed_report`: Generate a detailed report of the BLAST search results.
-   `visualize_hits`: Generate a graphical visualization of the hits on the query sequence.
-   `export_results`: Export the BLAST search results in the specified format.
-   `batch_blast_search`: Perform a batch BLAST search on a sequence file in FASTA format.
-   `custom_blast_search`: Perform a custom BLAST search with custom parameters.

### Use example: blast_search

The `blast_search` function performs a BLAST search on the input sequence using the selected program and database, and returns the results.

Here is an example of how to use the `blast_search` function:

    from blast_query_tool import blast_search
    
    sequence = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
    program = "blastn"
    database = "GenBank"
    evalue_cutoff = 10
    max_hits = 100
    filter = "none"
    
    blast_record = blast_search(sequence, program, database, evalue_cutoff, max_hits, filter)

For more information on the usage of each function, refer to the code documentation.

## Contributing

If you would like to contribute to the development of this tool, please feel free to fork the repository and submit a pull request.

## License

This project is licensed under the MIT License.
