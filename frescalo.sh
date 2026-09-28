#!/usr/bin/env bash
# Development/maintenance helper for Frescalo_rs.
# Run `./frescalo.sh help` for the list of subcommands.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

usage() {
  cat <<'EOF'
Usage: ./frescalo.sh <command> [args]

Commands:
  build             cargo build --release (add --debug for a debug build)
  test              cargo test --release
  doc [--open]      cargo doc --no-deps (optionally open it in a browser)
  fmt               cargo fmt
  lint              cargo clippy --all-targets -- -D warnings
  verify            run verify_against_fortran.sh (needs Original_frescalo/)
  check             fmt --check, lint, test (what CI runs), plus verify if
                    Original_frescalo/ is present locally
  clean             cargo clean, plus local verify_work/ and target/doc/
  release <version> bump Cargo.toml to <version>, commit, tag vX.Y.Z, push
  help              show this message
EOF
}

cmd_build() {
  if [ "${1:-}" = "--debug" ]; then
    cargo build
  else
    cargo build --release
  fi
}

cmd_test() {
  cargo test --release
}

cmd_doc() {
  cargo doc --no-deps
  if [ "${1:-}" = "--open" ]; then
    cargo doc --no-deps --open
  fi
}

cmd_fmt() {
  cargo fmt "$@"
}

cmd_lint() {
  cargo clippy --all-targets -- -D warnings
}

cmd_verify() {
  if [ ! -d Original_frescalo ]; then
    echo "error: Original_frescalo/ not found (it's gitignored; keep a local" >&2
    echo "copy of the original Fortran sources/data there to run this check)" >&2
    exit 1
  fi
  ./verify_against_fortran.sh
}

cmd_check() {
  cargo fmt --check
  cmd_lint
  cmd_test
  if [ -d Original_frescalo ]; then
    cmd_verify
  else
    echo "note: Original_frescalo/ not found locally; skipping verify (not part of CI)"
  fi
}

cmd_clean() {
  cargo clean
  rm -rf verify_work target/doc
}

cmd_release() {
  local version="${1:-}"
  if [ -z "$version" ]; then
    echo "usage: ./frescalo.sh release <version>   (e.g. 1.1.0)" >&2
    exit 1
  fi
  if [ -n "$(git status --porcelain)" ]; then
    echo "error: working tree is not clean; commit or stash first" >&2
    exit 1
  fi
  sed -i.bak "s/^version = \".*\"/version = \"${version}\"/" Cargo.toml
  rm -f Cargo.toml.bak
  cargo check --release >/dev/null # refresh Cargo.lock's version entry
  git add Cargo.toml Cargo.lock
  git commit -m "Release v${version}"
  git tag "v${version}"
  echo "Created commit and tag v${version}. Push with:"
  echo "  git push && git push origin v${version}"
}

case "${1:-help}" in
  build) shift; cmd_build "$@" ;;
  test) shift; cmd_test "$@" ;;
  doc) shift; cmd_doc "$@" ;;
  fmt) shift; cmd_fmt "$@" ;;
  lint) shift; cmd_lint "$@" ;;
  verify) shift; cmd_verify "$@" ;;
  check) shift; cmd_check "$@" ;;
  clean) shift; cmd_clean "$@" ;;
  release) shift; cmd_release "$@" ;;
  help|-h|--help) usage ;;
  *) echo "error: unknown command '${1}'" >&2; usage; exit 1 ;;
esac
