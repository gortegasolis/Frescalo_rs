//! Rust port of neighsim_1.f:
//! NEIGHSIM - Neighbourhood similarity based on training-set species and
//! physical proximity
//! written by Mark Hill, January-June 2011
//!
//! Calculates floristic similarity between samples from a training set,
//! combined with physical proximity ranks, and writes neighbourhood weights
//! for use in FRESCALO.

use frescalo::cli::{Args, OptSpec, Session, NO_HOLD};
use frescalo::config::{limits, NeighsimParams, DEFAULT_NEIGHBOURHOOD_SIZE};
use frescalo::fortran::{cout, ifmt, ld_i, ld_line, DataReader};
use frescalo::neighbourhood::similarity::{read_inputs, shared_species};
use frescalo::neighbourhood::weights::{mark_spatial_shortlist, write_neighbourhoods};
use std::io::{BufWriter, Write};

const PROGRAM: &str = "neighsim";

const OPTIONS: &[OptSpec] = &[
    OptSpec {
        name: "training",
        value: Some("FILE"),
        help: "Training-set species data [site species]",
    },
    OptSpec {
        name: "distances",
        value: Some("FILE"),
        help: "Neighbourhood distances from sampdist",
    },
    OptSpec {
        name: "similarity-out",
        value: Some("FILE"),
        help: "Training-set similarity output file (must not exist)",
    },
    OptSpec {
        name: "weights-out",
        value: Some("FILE"),
        help: "Neighbourhood weights output file for frescalo (must not exist)",
    },
    OptSpec {
        name: "neighbours",
        value: Some("K"),
        help: "Most similar sites kept per site [default: 100]",
    },
    NO_HOLD,
];

fn main() {
    let args = Args::parse(
        PROGRAM,
        "NEIGHSIM - Neighbourhood similarity based on training-set species and physical proximity",
        OPTIONS,
    );
    let mut session = Session::new(PROGRAM, args.has("no-hold"));

    cout("");
    cout(" NEIGHSIM - Neighbourhood similarity based on training-set species and physical proximity");
    cout(" written by Mark Hill, January-June 2011");
    cout("");
    // 1x,'NOTE PARAMETER LIMITS: Sites',i5,' Species',i6,1x,'Training-set number of records',i9
    cout(&format!(
        " NOTE PARAMETER LIMITS: Sites{} Species{} Training-set number of records{}",
        ifmt(limits::NEIGHSIM_SITES as i64, 5),
        ifmt(limits::NEIGHSIM_SPECIES as i64, 6),
        ifmt(limits::NEIGHSIM_RECORDS as i64, 9)
    ));
    cout("");
    let (_f1, ftr) = session.input_file(
        args.get("training"),
        &[" Type name of input file with Training-set species data [sample species] ...."],
    );
    let mut train_reader = DataReader::new(ftr);
    let (_f2, fdi) = session.input_file(
        args.get("distances"),
        &[" Type name of input file with physical distances ...."],
    );
    let mut dist_reader = DataReader::new(fdi);
    let (_f3, fsim) = session.output_file(
        args.get("similarity-out"),
        &[" Type name of Training-set similarity output file ..."],
    );
    let mut unit9 = BufWriter::new(fsim);
    let (_f4, fwgt) = session.output_file(
        args.get("weights-out"),
        &[" Type name of weights output file for use in Frescalo ..."],
    );
    let mut unit8 = BufWriter::new(fwgt);
    let params = NeighsimParams {
        neighbours: session.neighbours(&args, DEFAULT_NEIGHBOURHOOD_SIZE),
    };

    let inp = match read_inputs(&mut dist_reader, &mut train_reader) {
        Ok(inp) => inp,
        Err(msg) => {
            ld_line(&msg);
            session.fail();
        }
    };
    let mut simil = match shared_species(&inp) {
        Ok(s) => s,
        Err(msg) => {
            ld_line(&msg);
            session.fail();
        }
    };

    // Multiply those within the spatial shortlist by BIG
    let (iseqq, neigh1) = mark_spatial_shortlist(&inp, &mut simil);
    if neigh1 == 0 {
        // list-directed write to the weights file: leading blank
        unit8.write_all(b" Unrecognized sample names in distance data\n").unwrap();
    }

    write_neighbourhoods(&inp, &simil, &iseqq, params.neighbours, neigh1, &mut unit9, &mut unit8);
    ld_line(&format!(
        "neigh,neigh1{}{}",
        ld_i(params.neighbours as i64),
        ld_i(neigh1 as i64)
    ));

    unit8.flush().unwrap();
    unit9.flush().unwrap();
    session.finish();
}
