use assert_cmd::Command;
use predicates::prelude::*;
use std::io::Read;

const PANEL: &str = "tests/cases/panel.tsv";
const ANNOTATION: &str = "tests/cases/ann.gff3";
const REF: &str = "tests/cases/ref.fa";

#[test]
fn build_invalid_panel_file() -> Result<(), Box<dyn std::error::Error>> {
    let sub = "build";
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(&[sub, "-i file/doesnt/exist", "-a", "gff", "-f", "fa"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file"));

    Ok(())
}

#[test]
fn build_full_run() -> Result<(), Box<dyn std::error::Error>> {
    let sub = "build";
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    let outdir = tempfile::tempdir()?;
    cmd.args(&[
        sub,
        "-i",
        PANEL,
        "-a",
        ANNOTATION,
        "-f",
        REF,
        "-o",
        outdir.path().to_string_lossy().as_ref(),
    ]);
    cmd.assert().success();

    let mut file1 = std::fs::File::open("tests/cases/expected/dr.prg")?;
    let mut file2 =
        std::fs::File::open(format!("{}/dr.prg", outdir.path().to_string_lossy()))?;

    let mut contents = String::new();
    file1.read_to_string(&mut contents)?;
    let mut other = String::new();
    file2.read_to_string(&mut other)?;

    let mut sorted1 = contents.as_bytes().to_owned();
    sorted1.sort_unstable();

    let mut sorted2 = other.as_bytes().to_owned();
    sorted2.sort_unstable();

    assert_eq!(sorted1, sorted2);

    Ok(())
}
