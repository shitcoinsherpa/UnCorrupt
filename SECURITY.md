# Security Policy

## Supported versions

Only the latest minor release of UnCorrupt receives security fixes. Older releases are archived snapshots and will not be patched.

| Version | Status |
|---|---|
| 1.x | Supported |
| < 1.0 | Pre-release, unsupported |

## Reporting a vulnerability

Please **do not file public GitHub issues for security problems**. Posting a working exploit before a fix is available puts every user of the tool at risk.

The fastest, safest path is a private GitHub Security Advisory:

1. Go to https://github.com/shitcoinsherpa/UnCorrupt/security/advisories/new
2. Describe the issue. Include:
   - The smallest input file or command that demonstrates the problem
   - The version of `uncorrupt` you tested against (`uncorrupt --version`)
   - The Python and OS you ran it on
   - What you expected, and what actually happened
3. We will acknowledge within 72 hours, triage within one week, and target a fix release within 30 days for any confirmed vulnerability.

If GitHub Security Advisories are not available to you, open an issue titled "Security: please contact me" with no details, and we will reach out via the contact on your profile.

## What counts as a vulnerability here

UnCorrupt reads spreadsheets supplied by the user. The interesting attack surfaces are:

- **Spreadsheet parsing**: a crafted `.xlsx` / `.xls` / `.csv` that causes arbitrary file write, code execution, or sustained denial of service when passed to `uncorrupt detect`.
- **Schema export**: a crafted input that causes `uncorrupt schema` to emit a JSON sidecar that, when validated by Frictionless, leaks paths or hangs.
- **The Gradio UI**: any path traversal, file-disclosure, or persistent server state introduced by the upload/download flow.
- **Supply chain**: a published wheel, container image, or release artifact whose contents do not match the source at the corresponding git tag.

Bug reports about misdetected gene symbols are **not** security issues. Please file those as normal GitHub issues with the offending fixture attached.

## Verification

Releases are reproducible from the source at the matching tag. You can verify:

- **Wheel**: `pip download --no-deps uncorrupt==<version>`, then `python -m build` from the corresponding tag and compare the resulting wheel's SHA256.
- **Container**: `docker pull ghcr.io/shitcoinsherpa/uncorrupt:<version>` and inspect the OCI labels for `org.opencontainers.image.source` and `org.opencontainers.image.version`.
- **SBOM**: each GitHub Release attaches a CycloneDX SBOM (`sbom.json`).
