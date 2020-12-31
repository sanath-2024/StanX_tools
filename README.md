# Stan-X Tools
Various tools associated with mapping transposons (originated as part of the StanEx project)

As of December 2020, a transposable element mapper, variant caller, synthetic reference genome generator, and download helper are provided.

These tools are subcommands of the main "sx" command. For example:

```bash
sx map [ARGS]
sx variants [ARGS]
sx sg [ARGS]
sx download [ARGS]
```

Use the `--help` argument for more information on the different subcommands and arguments:

```bash
sx --help
```

Note: the TE mapper uses an algorithm that is based on Bergman Lab's `ngs_te_mapper` tool, written by Raquel S. Linheiro, Michael G. Nelson, and Casey M. Bergman.

- ngs_te_mapper link: https://github.com/bergmanlab/ngs_te_mapper
- ngs_te_mapper license: https://github.com/bergmanlab/ngs_te_mapper/blob/master/LICENSE