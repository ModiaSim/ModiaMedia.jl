using Documenter, ModiaMedia

makedocs(
  modules  = [ModiaMedia],
  format   = :html,
  sitename = "ModiaMedia",
  authors  = "Martin Otter (DLR-SR), Hilding Elmqivst (Mogram), Chris Laughman (NERL)",
  html_prettyurls = false,
  pages    = [
     "Home"   => "index.md",
     "Library" => [
        "lib/MediaTypes.md",
        "lib/Types.md",
        "lib/Functions.md"
        ]
  ]
)

