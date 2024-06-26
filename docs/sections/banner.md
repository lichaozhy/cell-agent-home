<v-responsive
  :aspect-ratio="16 / 9"
  class="border-0 px-md-16 px-sm-0 py-0"
>
  <template v-if="isCN !== null">
    <iframe
      v-if="isCN"
      src="//player.bilibili.com/player.html?isOutside=true&aid=112681956672810&bvid=BV14Z3MexEoy&cid=500001596752173&p=1&rel=0"
      scrolling="no"
      allowfullscreen="true"
      class="h-100 w-100 border-0"
    ></iframe>
    <iframe
      v-else
      class="h-100 w-100 border-0"
      src="https://www.youtube.com/embed/kJONeY-DDeE?si=Gl0P4QR2Avs9CocD&rel=0"
      title="YouTube video player"
      allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share"
      referrerpolicy="strict-origin-when-cross-origin"
      allowfullscreen
      ></iframe>
  </template>
</v-responsive>
