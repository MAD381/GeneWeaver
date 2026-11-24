# GeneWeaver Project – Weekly Progress Report

## Week 1: FASTA Processing and DNA Sequence Analysis

### Objective
The main objective of Week 1 was to develop the foundation of the GeneWeaver project by implementing a FASTA processing pipeline capable of reading, validating, and analyzing DNA sequences.

### Work Completed

1. **FASTA File Parsing**
   - Implemented FASTA file reading using the Biopython library.
   - Extracted DNA sequences and their corresponding identifiers.
   - Added handling for missing or empty FASTA files.

2. **DNA Sequence Validation**
   - Implemented validation to identify invalid DNA characters.
   - Added support for lowercase DNA sequences.
   - Handled empty sequences and invalid input conditions.

3. **DNA Sequence Chunking**
   - Implemented a sequence chunking mechanism to divide large DNA sequences into smaller sections.
   - Added configurable chunk size and overlap parameters.
   - Ensured overlapping regions are maintained between consecutive chunks.

4. **Sequence Analysis**
   - Implemented basic DNA sequence statistics.
   - Added GC-content calculation and sequence fingerprinting.
   - Implemented reverse complement generation.
   - Added duplicate sequence detection.

### Outcome
Successfully developed the core FASTA processing functionality, providing a foundation for further DNA analysis and processing in GeneWeaver.

---

## Week 2: CLI Integration and Processing Improvements

### Objective
The objective of Week 2 was to improve the FASTA processing workflow by integrating it with a command-line interface (CLI), adding progress monitoring, processing-time measurement, and better input validation.

### Work Completed

1. **Command-Line Interface Integration**
   - Integrated the FASTA processing module with the project's `cli.py`.
   - Added command-line arguments for input file, output directory, chunk size, and overlap.
   - Enabled users to execute the pipeline without manually modifying the source code.

2. **Input Validation and Error Handling**
   - Added validation for positive chunk sizes.
   - Added validation to prevent negative overlap values.
   - Ensured overlap is smaller than the chunk size.
   - Added handling for missing input files and processing errors.

3. **Progress Monitoring**
   - Integrated the TQDM progress bar.
   - Displayed the number of sequences processed.
   - Displayed chunk-generation progress during execution.

4. **Processing-Time Measurement**
   - Integrated execution-time measurement using Python's `time.perf_counter()`.
   - Displayed the total processing time after successful execution.

5. **Output Generation**
   - Generated processed FASTA output files.
   - Organized generated files inside the specified output directory.

6. **Automated Testing**
   - Added automated CLI tests in `tests/test_cli.py`.
   - Tested valid inputs, invalid parameters, missing files, successful processing, and progress callback integration.
   - Executed the complete unit test suite.

### Outcome
Successfully integrated the FASTA processing module with the CLI and improved usability through progress monitoring, error handling, and processing-time reporting.

---

## Technologies and Tools Used

- Python
- Biopython
- TQDM
- Argparse
- Unittest
- Git and GitHub
- Visual Studio Code

## Testing Results

| Testing Category | Result |
|---|---:|
| Existing project tests | 49 passed |
| New CLI tests | 8 passed |
| Total automated tests | 57 passed |
| CLI execution | Successful |
| FASTA output generation | Successful |

## Overall Conclusion

During Weeks 1 and 2, the core FASTA processing pipeline was developed and integrated with the GeneWeaver command-line interface. The implementation supports DNA sequence parsing, validation, chunking, basic sequence analysis, progress monitoring, and output generation.

The project currently has 57 passing automated tests, demonstrating that the implemented functionality works correctly under the tested conditions. This work establishes a foundation for the next stages of GeneWeaver development and integration.