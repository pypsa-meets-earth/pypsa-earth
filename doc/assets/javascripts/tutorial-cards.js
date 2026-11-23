// SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
//
// SPDX-License-Identifier: AGPL-3.0-or-later

// Make each tutorial landing page card clickable: clicking anywhere on the
// panel (outside of a specific link) opens the first tutorial page in that card.
document.addEventListener("DOMContentLoaded", function () {
  var cards = document.querySelectorAll(".tutorial-cards > ul > li");

  cards.forEach(function (card) {
    var firstLink = Array.prototype.find.call(
      card.querySelectorAll("a"),
      function (link) { return !link.closest(".tutorial-image-source"); }
    );
    if (!firstLink) return;

    card.addEventListener("click", function (event) {
      if (event.target.closest("a")) return;
      window.location.href = firstLink.getAttribute("href");
    });
  });
});
