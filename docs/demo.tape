# UnCorrupt — README demo
#
# Run:  vhs docs/demo.tape   (https://github.com/charmbracelet/vhs)
# Output: docs/demo.gif

Output docs/demo.gif

Set FontSize 16
Set Width 1200
Set Height 700
Set Theme "Catppuccin Mocha"
Set TypingSpeed 50ms

Type "# Scan a supplementary spreadsheet for Excel-style gene-symbol corruption"
Enter
Sleep 800ms

Type "uncorrupt detect tests/fixtures/fixture_01_day_of_month.xlsx"
Enter
Sleep 3s

Type "clear"
Enter
Sleep 300ms

Type "# Emit a Frictionless schema sidecar that prevents re-corruption on import"
Enter
Sleep 800ms

Type "uncorrupt schema tests/fixtures/fixture_01_day_of_month.xlsx"
Enter
Sleep 2500ms

Type "clear"
Enter
Sleep 300ms

Type "# Verify the schema rejects a corrupted copy"
Enter
Sleep 800ms

Type "frictionless validate --schema fixture_01_day_of_month.schema.json corrupted_copy.csv"
Enter
Sleep 2500ms

Type "# 100% post-boost precision (Wilson 95% [0.998, 1.000])"
Enter
Sleep 800ms
Type "# 98.21% recall on Ziemann S2"
Enter
Sleep 800ms
Type "# 100% recall on independent Ziemann 2016 corpus"
Enter
Sleep 2s
