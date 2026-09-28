//! Command-line options, with the original interactive prompts as fallback.
//!
//! Rules, shared by all three programs:
//!
//! * **No options:** the program is fully interactive and asks every question
//!   the Fortran original asks, in the same order. A blank answer to a
//!   question that has a default (number of neighbours, Φ, R\*) takes the
//!   default.
//! * **Any value option given:** each setting comes from its option. A file
//!   name that is not given is still asked for at the prompt; settings with a
//!   default (neighbours, Φ, R\*, exclusion file, distance method) just take
//!   the default without asking.
//! * The final "Press \<RETURN\> to exit" pause happens only if something was
//!   read from the keyboard; `--no-hold` suppresses it in every case.
//!
//! The parser is deliberately small (`--name value` or `--name=value`, and
//! switches), which keeps the project free of dependencies that would need a
//! newer compiler than the one the port is built with.

use crate::fortran::console::{cout, filin, filout, hold, read_a20, read_line};
use std::fs::{File, OpenOptions};
use std::io::{self, StdinLock, Write};
use std::str::FromStr;

/// One command-line option.
pub struct OptSpec {
    /// Long name, without the leading `--`.
    pub name: &'static str,
    /// Placeholder for the value (`None` for a switch).
    pub value: Option<&'static str>,
    /// One-line help text.
    pub help: &'static str,
}

/// `--no-hold`, accepted by every program.
pub const NO_HOLD: OptSpec = OptSpec {
    name: "no-hold",
    value: None,
    help: "Do not wait for <RETURN> before exiting",
};

/// Parsed command line.
pub struct Args {
    values: Vec<(&'static str, String)>,
    switches: Vec<&'static str>,
}

impl Args {
    /// Parse `std::env::args()` against `specs`. `--help` prints usage and
    /// exits; unknown or malformed options print an error and exit with
    /// status 2.
    pub fn parse(program: &str, about: &str, specs: &[OptSpec]) -> Args {
        Self::parse_from(program, about, specs, std::env::args().skip(1).collect())
    }

    fn parse_from(program: &str, about: &str, specs: &[OptSpec], argv: Vec<String>) -> Args {
        let mut args = Args {
            values: Vec::new(),
            switches: Vec::new(),
        };
        let mut it = argv.into_iter();
        while let Some(a) = it.next() {
            if a == "--help" || a == "-h" {
                print_usage(program, about, specs);
                std::process::exit(0);
            }
            let body = match a.strip_prefix("--") {
                Some(b) => b,
                None => usage_error(program, &format!("unexpected argument '{}'", a)),
            };
            let (name, inline) = match body.split_once('=') {
                Some((n, v)) => (n, Some(v.to_string())),
                None => (body, None),
            };
            let spec = match specs.iter().find(|s| s.name == name) {
                Some(s) => s,
                None => usage_error(program, &format!("unknown option '--{}'", name)),
            };
            if spec.value.is_some() {
                let v = match inline.or_else(|| it.next()) {
                    Some(v) => v,
                    None => usage_error(program, &format!("option '--{}' needs a value", name)),
                };
                args.values.retain(|(n, _)| *n != spec.name);
                args.values.push((spec.name, v));
            } else {
                if inline.is_some() {
                    usage_error(program, &format!("option '--{}' takes no value", name));
                }
                args.switches.push(spec.name);
            }
        }
        args
    }

    /// Value of option `name`, if given.
    pub fn get(&self, name: &str) -> Option<&str> {
        self.values.iter().find(|(n, _)| *n == name).map(|(_, v)| v.as_str())
    }

    /// Whether switch `name` was given.
    pub fn has(&self, name: &str) -> bool {
        self.switches.iter().any(|n| *n == name)
    }

    /// True if any value option was given (see the module rules).
    pub fn option_mode(&self) -> bool {
        !self.values.is_empty()
    }

    /// Value of option `name` parsed as `T`; exits with status 2 if it does
    /// not parse.
    pub fn parse_value<T: FromStr>(&self, program: &str, name: &str) -> Option<T> {
        self.get(name).map(|v| match v.parse::<T>() {
            Ok(x) => x,
            Err(_) => usage_error(program, &format!("invalid value '{}' for '--{}'", v, name)),
        })
    }
}

fn print_usage(program: &str, about: &str, specs: &[OptSpec]) {
    println!("{}\n\nUsage: {} [OPTIONS]\n\nOptions:", about, program);
    for s in specs {
        let left = match s.value {
            Some(v) => format!("--{} <{}>", s.name, v),
            None => format!("--{}", s.name),
        };
        println!("  {:<30} {}", left, s.help);
    }
    println!("  {:<30} {}", "-h, --help", "Print this help");
    println!("\nWith no options the program asks for everything interactively.");
}

fn usage_error(program: &str, msg: &str) -> ! {
    eprintln!("{}: {}\nTry '{} --help'.", program, msg, program);
    std::process::exit(2);
}

/// Print an error for a run driven by options and exit with status 1.
pub fn fatal(program: &str, msg: &str) -> ! {
    eprintln!("{}: {}", program, msg);
    std::process::exit(1);
}

/// Console session: reads prompted answers from stdin and remembers whether
/// anything was read, which decides whether the final pause happens.
pub struct Session {
    pub program: &'static str,
    stdin: StdinLock<'static>,
    prompted: bool,
    no_hold: bool,
}

impl Session {
    pub fn new(program: &'static str, no_hold: bool) -> Self {
        Session {
            program,
            stdin: io::stdin().lock(),
            prompted: false,
            no_hold,
        }
    }

    /// Raw access to stdin for the Fortran-style readers, marking the session
    /// as interactive.
    pub fn stdin(&mut self) -> &mut StdinLock<'static> {
        self.prompted = true;
        &mut self.stdin
    }

    /// Open an existing input file: the given name, or ask at the prompt
    /// (after printing `prompt`) until an existing file is named.
    pub fn input_file(&mut self, given: Option<&str>, prompt: &[&str]) -> (String, File) {
        match given {
            Some(name) => match File::open(name) {
                Ok(f) => (name.to_string(), f),
                Err(_) => fatal(self.program, &format!("input file '{}' does not exist", name)),
            },
            None => {
                for p in prompt {
                    cout(p);
                }
                filin(self.stdin())
            }
        }
    }

    /// Create a new output file (never overwriting): the given name, or ask
    /// at the prompt.
    pub fn output_file(&mut self, given: Option<&str>, prompt: &[&str]) -> (String, File) {
        match given {
            Some(name) => match OpenOptions::new().write(true).create_new(true).open(name) {
                Ok(f) => (name.to_string(), f),
                Err(_) => fatal(
                    self.program,
                    &format!("output file '{}' already exists or cannot be created", name),
                ),
            },
            None => {
                for p in prompt {
                    cout(p);
                }
                filout(self.stdin())
            }
        }
    }

    /// Read a file name (first 20 characters) at the prompt.
    pub fn read_a20(&mut self) -> String {
        read_a20(self.stdin())
    }

    /// Number of neighbours: the option value, the default in option mode,
    /// or the prompt (blank answer = default). Must be at least 1.
    pub fn neighbours(&mut self, args: &Args, default: i32) -> i32 {
        if let Some(v) = args.parse_value::<i32>(self.program, "neighbours") {
            if v < 1 {
                fatal(self.program, "--neighbours must be at least 1");
            }
            return v;
        }
        if args.option_mode() {
            return default;
        }
        loop {
            cout(" Type number of neighbours to include ...");
            let line = read_line(self.stdin());
            let tok = match line.split_whitespace().next() {
                None => return default,
                Some(t) => t,
            };
            match tok.parse::<i32>() {
                Ok(v) if v >= 1 => return v,
                Ok(_) => cout(" ***ERROR*** Must be at least 1"),
                Err(_) => {
                    eprintln!("\nFortran runtime error: Bad integer for item 1 in list input");
                    std::process::exit(2);
                }
            }
        }
    }

    /// Normal end of run: pause for `<RETURN>` if interactive, then exit 0.
    pub fn finish(&mut self) -> ! {
        let _ = io::stdout().flush();
        if self.prompted && !self.no_hold {
            hold(&mut self.stdin);
        }
        std::process::exit(0);
    }

    /// Fatal data error (message already printed): pause if interactive, as
    /// the originals do, then exit with status 1.
    pub fn fail(&mut self) -> ! {
        let _ = io::stdout().flush();
        if self.prompted && !self.no_hold {
            hold_then(&mut self.stdin, 1);
        }
        std::process::exit(1);
    }
}

fn hold_then(stdin: &mut StdinLock<'static>, code: i32) -> ! {
    let mut o = io::stdout();
    let _ = o.write_all(b"\n\nPress <RETURN> to exit\n\n\n----------------------\n");
    let _ = o.flush();
    let mut s = String::new();
    let _ = io::BufRead::read_line(stdin, &mut s);
    std::process::exit(code);
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPECS: &[OptSpec] = &[
        OptSpec { name: "neighbours", value: Some("N"), help: "" },
        OptSpec { name: "out", value: Some("FILE"), help: "" },
        NO_HOLD,
    ];

    fn parse(v: &[&str]) -> Args {
        Args::parse_from("t", "t", SPECS, v.iter().map(|s| s.to_string()).collect())
    }

    #[test]
    fn values_and_switches() {
        let a = parse(&["--neighbours", "50", "--out=x.txt", "--no-hold"]);
        assert_eq!(a.get("neighbours"), Some("50"));
        assert_eq!(a.get("out"), Some("x.txt"));
        assert!(a.has("no-hold"));
        assert!(a.option_mode());
        assert_eq!(a.parse_value::<i32>("t", "neighbours"), Some(50));
    }

    #[test]
    fn switch_alone_keeps_interactive_mode() {
        let a = parse(&["--no-hold"]);
        assert!(!a.option_mode());
        assert_eq!(a.get("neighbours"), None);
    }

    #[test]
    fn last_value_wins() {
        let a = parse(&["--neighbours", "5", "--neighbours", "7"]);
        assert_eq!(a.get("neighbours"), Some("7"));
    }
}
