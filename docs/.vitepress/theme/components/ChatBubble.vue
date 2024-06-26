<template>
  <v-sheet class="mt-10 mb-10 position-relative" border="sm" rounded>
    <v-avatar
      size="28"
      class="chat-bubble-subject position-absolute"
    >
      <slot name="subject">
      </slot>
    </v-avatar>

    <div class="pa-2">
      <slot>
        <div>content</div>
      </slot>
    </div>
  </v-sheet>

</template>

<script setup>
import { mdiAccount } from '@mdi/js'
</script>

<style>
.chat-bubble-subject {
  top: -32px;
}
</style>
