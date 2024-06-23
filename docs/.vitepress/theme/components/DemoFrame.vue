<template>
  <v-container class="px-0">
    <v-row>
      <v-col cols="6">
        <slot name="description">Description</slot>
      </v-col>
      <v-col cols="6">
        <v-card
          border="surface-variant sm opacity-100"
          title="Input"
          variant="text"
          class="rounded-0"
        >
          <v-card-text>
            <slot name="input">Input</slot>
          </v-card-text>
        </v-card>
        <v-card
          border="surface-variant sm opacity-100"
          title="Thought"
          variant="text"
          class="rounded-0"
        >
          <v-card-text>
            <slot name="thought">Thought</slot>
          </v-card-text>
        </v-card>
        <v-card
          border="surface-variant sm opacity-100"
          title="Output"
          variant="text"
          class="rounded-0 border-t-0"
        >
          <v-card-text>
            <slot name="output">Output</slot>
					</v-card-text> </v-card
      ></v-col>
    </v-row>
  </v-container>
</template>

<script setup></script>
