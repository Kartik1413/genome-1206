document.addEventListener('DOMContentLoaded', function() {
    const moreInfoLinks = document.querySelectorAll('.scroll-to-section a');

    moreInfoLinks.forEach(link => {
        link.addEventListener('click', function(event) {
            event.preventDefault();
            window.location.href = this.getAttribute('href');
        });
    });
});
