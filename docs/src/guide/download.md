# Download prebuilt index

The `index` subcommand can be used to download or list available indices.

To download the latest index version for all species

```
drprg index --download
# is the same as
drprg index --download all
```

If you only want the latest for a particular species

```
drprg index --download mtb
# is the same as
drprg index --download mtb@latest
```

If you want to see what indices are available in the first place

```
drprg index --list
+--------------+---------+----------+------------+
| Name         | Species | Version  | Downloaded |
+--------------+---------+----------+------------+
| mtb@20230308 | mtb     | 20230308 | N          |
+--------------+---------+----------+------------+
```

You can specify a version with the following syntax

```
drprg index --download mtb@<version>
```

where `<version>` is the version you want. If no version is provided, `latest` is used.

## Output directory

By default, indices are stored in `$HOME/.drprg/`, nested as `species/species-version`.
So downloading version `20230308` of species `mtb` will produce an index
at `$HOME/.drprg/mtb/mtb-20230308/`.

You can change the default output directory (`$HOME/.drprg/`) with the `--outdir`
option. If you do this, you'll need to pass the full path to the index when
using [`predict`](predict.md).

## Full usage

```
$ drprg index --help
Download and interact with indices

Usage: drprg index [OPTIONS] [NAME]

Arguments:
  [NAME]
          The name/path of the index to download

          [default: all]

Options:
  -d, --download
          Download a prebuilt index

  -v, --verbose
          Use verbose output

  -l, --list
          List all available (and downloaded) indices

  -t, --threads <INT>
          Maximum number of threads to use

          Use 0 to select the number automatically

          [default: 1]

  -o, --outdir <DIR>
          Index directory

          Use this if your indices are not in a default location, or you want to download them to a non-default location

          [default: $HOME/.drprg/]

  -F, --force
          Overwrite any existing indices

  -h, --help
          Print help (see a summary with '-h')
```