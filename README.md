# Cryptosporidium GP60 Characterization Tool


### Tool Functionality and Overview

This tool's intended purpose is to perform Cryptosporidium subtyping utilizing the GP60 genomic target region within the Cryptosporidum genome. The region targets the subtype family, short tandem repeats and secondary repeats within the sequece. The tool works with Sanger sequences and isolate whole genome contig sequences in fasta format. 

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

**Tool disclaimer**  ***Please note that the assays used are not ISO or CLIA-certified and should NOT be considered diagnostic!***


## Building Container

To rebuild the container, clone this repository and run the command below:
`docker build -tag(optional) -file Dockerfile <location>`

## Requirements:
- Linux Operating system 
- [Docker](https://docs.docker.com/) 
- Perl

## Expected Input:
- --fasta		Fasta file with 1 or more sanger sequences OR a fasta file containing the assembly of 1 samples
- --blastdb	Prebuilt curated blast database
- --data		Data type can be sanger or wgs (capitalization does not matter)

## Expected Output: 
Generates the gp60 subtyping results in a tab delimited text format

### Running with Docker
``` docker run -v $(pwd)/SM_TestData:/test --privileged wdpbcdsphl/crypto_gp60:2.4 perl /scripts/gp60Typer.pl --blastdb db/Crypto_GP60_DB --fasta /test/testInput_SM.fasta --data sanger > results_gp60SM-15.txt ```


### Running with Singularity
```singularity exec -B $(pwd)/SM_TestData:/dataIn crypto_gp60-2.4.simg perl /scripts/gp60Typer.pl --blastdb /db/Crypto_GP60_DB --fasta /dataIn/gp60_seqs.fasta --data sanger```

## Developer
Developed by: Alyssa Kelly, akelley139@gmail.com

Modification by: Shatavia Morrison, SMorrison@cdc.gov

Manuscript prepared by: Anusha Ginni, qxu0@cdc.gov

Clinical Detection Surveillance/WDPB, CDC

Database was updated on: 2021-10-26








## Related documents

* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)



  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
