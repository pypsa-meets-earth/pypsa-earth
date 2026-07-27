// Make each tutorial landing page card clickable: clicking anywhere on the
// panel (outside of a specific link) opens the first tutorial page in that card.
document.addEventListener("DOMContentLoaded", function () {
  var cards = document.querySelectorAll(".tutorial-cards > ul > li");

  cards.forEach(function (card) {
    var firstLink = card.querySelector("a");
    if (!firstLink) return;

    card.addEventListener("click", function (event) {
      if (event.target.closest("a")) return;
      window.location.href = firstLink.getAttribute("href");
    });
  });
});
