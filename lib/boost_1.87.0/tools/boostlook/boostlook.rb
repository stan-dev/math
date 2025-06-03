Asciidoctor::Extensions.register do
  postprocessor do
    process do |doc, output|
      output = output.sub(/(<body[^>]*>)/, '\1<div class="boostlook">')
      output = output.sub('</body>', '</div></body>')
      output = output.sub(/(<body.*?<div[^>]*id="toc"[^>]*>)/m, '\1<button id="toggle-toc" title="Show Table of Contents" aria-expanded="false" aria-controls="toc">☰</button>')
      output = output.sub(/(<body.*?<div[^>]*id="footer"[^>]*>)/m, '</div>\1')

      script_tag = <<~SCRIPT
        <script>
        (function() {
          const html = document.documentElement;
          const isPinned = localStorage.getItem('tocPinned') === 'true';

          html.classList.add('toc-hidden');
          if (isPinned) {
            html.classList.add('toc-pinned');
            html.classList.add('toc-visible');
            html.classList.remove('toc-hidden');
          }

          document.addEventListener("DOMContentLoaded", () => {
            const tocButton = document.getElementById("toggle-toc");
            const toc = document.getElementById("toc");

            if (!tocButton || !toc) return;

            let isPinned = localStorage.getItem('tocPinned') === 'true';

            function updateTocVisibility(visible) {
              html.classList.toggle("toc-visible", visible);
              html.classList.toggle("toc-hidden", !visible);
              tocButton.setAttribute("aria-expanded", visible);
              tocButton.textContent = visible ? "×" : "☰";
              tocButton.setAttribute("title", visible ? "Hide Table of Contents" : "Show Table of Contents");
            }

            tocButton.addEventListener("click", () => {
              isPinned = !isPinned;
              localStorage.setItem('tocPinned', isPinned);
              html.classList.toggle('toc-pinned', isPinned);
              updateTocVisibility(isPinned);  state
            });

            tocButton.addEventListener("mouseenter", () => {
              if (!isPinned) {
                updateTocVisibility(true);
              }
            });

            toc.addEventListener("mouseleave", () => {
              if (!isPinned) {
                updateTocVisibility(false);
              }
            });

            updateTocVisibility(isPinned);
          });
        })();
        </script>
      SCRIPT

      output = output.sub('</body>', "#{script_tag}</body>")

      output
    end
  end
end
