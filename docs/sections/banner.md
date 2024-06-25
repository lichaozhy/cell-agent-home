<v-responsive
  :aspect-ratio="16 / 9"
  class="border-0 px-md-16 px-sm-0 py-0"
>
  <iframe
    v-if="isCN"
    src="//player.bilibili.com/player.html?isOutside=true&aid=112613522411165&bvid=BV1dVGoeCEQ4&cid=500001581492325&p=1"
    scrolling="no"
    allowfullscreen="true"
    class="h-100 w-100 border-0"
  ></iframe>
  <iframe
    v-else
    class="h-100 w-100 border-0"
    src="https://www.youtube.com/embed/7a4M3ymp5ng?si=Hp-jAv9KkYHy-4-w&rel=0"
    title="YouTube video player"
    allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share"
    referrerpolicy="strict-origin-when-cross-origin"
    allowfullscreen
    ></iframe>
</v-responsive>
