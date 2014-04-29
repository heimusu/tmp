var getHotEntry;

getHotEntry = function() {
  return $.ajax({
    url: "http://feeds.feedburner.com/hatena/b/hotentry",
    dataType: "xml",
    success: function(data) {
      return $(data).find("item").each(function() {
        var title;
        title = $(this).find("title").text();
        return $(".feeds").append("<p>" + title + "</p>");
      });
    }
  });
};

$("button").on("click", function() {
  return getHotEntry();
});
